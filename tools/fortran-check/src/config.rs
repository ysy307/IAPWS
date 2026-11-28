use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet, HashMap},
    fs,
    path::{Path, PathBuf},
    sync::OnceLock,
};

use anyhow::{Context, Result, bail, ensure};
use globset::{GlobBuilder, GlobMatcher};
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    diagnostic::{Diagnostic, Severity, slash_path},
    library_calls::SafetyPolicy,
};

#[derive(Debug, Clone, Serialize)]
pub struct RuleSpec {
    pub id: &'static str,
    pub phase: u8,
    pub name: &'static str,
    pub default_severity: &'static str,
}

macro_rules! rules {
    ($($id:literal, $phase:literal, $name:literal, $severity:literal);+ $(;)?) => {
        &[$(RuleSpec { id: $id, phase: $phase, name: $name, default_severity: $severity }),+]
    };
}

pub const RULES: &[RuleSpec] = rules![
    "STD001", 1, "complete IMPLICIT NONE", "error";
    "STD002", 1, "INTENT(OUT)", "error";
    "STD003", 1, "explicit procedure interface", "error";
    "STD004", 1, "named kind", "error";
    "STD005", 1, "implicit SAVE", "error";
    "STD006", 1, "prohibited statement", "error";
    "C121", 1, "USE with an explicit ONLY list", "error";
    "OOP001", 2, "default module access", "error";
    "OOP002", 2, "public type component", "error";
    "OOP003", 2, "mutable module data", "error";
    "OOP004", 2, "explicit binding access", "error";
    "OOP005", 2, "public non-type-bound procedure", "error";
    "OOP006", 2, "exported binding implementation", "error";
    "OOP007", 2, "public module generic", "error";
    "OOP008", 2, "unlimited polymorphism", "error";
    "OOP009", 2, "NOPASS", "error";
    "OOP010", 2, "pointer data component", "error";
    "OOP011", 2, "generic specifics", "error";
    "SEM001", 3, "duplicate symbol", "error";
    "SEM002", 3, "unresolved symbol", "error";
    "SEM003", 3, "ambiguous or inaccessible symbol", "error";
    "SEM004", 3, "separate procedure or submodule", "error";
    "SEM005", 3, "unresolved binding or interface", "error";
    "DIT001", 4, "inheritance cycle", "error";
    "DIT002", 4, "unmeasurable base type", "error";
    "DIT003", 4, "inheritance depth", "error";
    "COH001", 4, "LCOM4", "error";
    "CBO001", 4, "coupling between types", "warning";
    "WMC001", 4, "weighted methods per type", "warning";
    "CC001", 4, "cyclomatic complexity", "warning";
    "ARCH001", 5, "source layer", "error";
    "ARCH002", 5, "duplicate module", "error";
    "ARCH003", 5, "unresolved module", "error";
    "ARCH004", 5, "forbidden layer edge", "error";
    "ARCH005", 5, "internal module access", "error";
    "ARCH006", 5, "external library access", "error";
    "ARCH007", 5, "module dependency cycle", "error";
    "MEM001", 6, "prefer allocatable", "warning";
    "MEM002", 6, "undefined pointer association", "conditional";
    "MEM003", 6, "null or released pointer storage", "conditional";
    "MEM004", 6, "implicit reallocation", "warning";
    "MEM005", 6, "array packing temporary", "warning";
    "OWN001", 6, "resource not released", "warning";
    "OWN002", 6, "resource release mismatch", "conditional";
    "OWN003", 6, "owning value copy", "warning";
    "OWN004", 6, "ownership cleanup or finalizer", "conditional";
    "OWN005", 6, "overwritten resource", "warning";
    "NUM001", 7, "numeric kernel purity", "conditional";
    "NUM002", 7, "external I/O in a numeric kernel", "error";
    "NUM003", 7, "mixed numeric kinds", "warning";
    "NUM004", 7, "implicit precision loss", "error";
    "NUM005", 7, "dynamic dispatch in a hot loop", "warning";
    "NUM006", 7, "allocation in a hot loop", "warning";
    "PRL001", 7, "MPI error status", "conditional";
    "PRL002", 7, "PETSc error status", "conditional";
    "PRL003", 7, "unfinished request or Begin/End", "warning";
    "PRL004", 7, "unsafe request operation", "conditional";
    "PRL005", 7, "rank-divergent collectives", "warning";
    "PRL006", 7, "parallel storage side effect", "conditional";
    "PAR001", 8, "unparseable Fortran", "error";
    "PAR002", 8, "unexpanded directive or include", "error";
];

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Limits {
    pub dit: usize,
    pub lcom4: usize,
    pub cbo: usize,
    pub wmc: usize,
    pub complexity: usize,
    pub subjects: BTreeMap<String, SubjectLimits>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubjectLimits {
    pub cbo: Option<usize>,
    pub lcom4: Option<usize>,
    pub wmc: Option<usize>,
    pub complexity: Option<usize>,
    pub reason: String,
}

impl Limits {
    pub fn complexity_for(&self, name: &str) -> usize {
        self.subjects
            .get(name)
            .and_then(|limits| limits.complexity)
            .unwrap_or(self.complexity)
    }
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            dit: 3,
            lcom4: 1,
            cbo: 5,
            wmc: 50,
            complexity: 10,
            subjects: BTreeMap::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Metrics {
    pub exclude_abstract: bool,
    pub exclude_stateless: bool,
    pub external_fan_out: bool,
    pub complexity_nodes: Vec<String>,
}

impl Default for Metrics {
    fn default() -> Self {
        Self {
            exclude_abstract: true,
            exclude_stateless: true,
            external_fan_out: true,
            complexity_nodes: list(&[
                "if_statement",
                "elseif_clause",
                "do_loop",
                "where_statement",
                "forall_statement",
                "case_statement",
                "type_statement",
                "rank_statement",
            ]),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Language {
    pub standard: String,
    pub line_length: usize,
    pub implicit_scopes: Vec<String>,
    pub require_implicit_external: bool,
    pub named_kinds: Vec<String>,
    pub banned_statements: Vec<String>,
}

impl Default for Language {
    fn default() -> Self {
        Self {
            standard: "f2018".into(),
            line_length: 132,
            implicit_scopes: list(&[
                "module",
                "submodule",
                "program",
                "function",
                "subroutine",
                "module_procedure",
            ]),
            require_implicit_external: true,
            named_kinds: list(&["integer", "real", "complex", "double"]),
            banned_statements: list(&["goto"]),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Numerical {
    pub kernels: Vec<String>,
    pub hot_loops: Vec<String>,
    pub purity_warnings: Vec<String>,
    pub io_statements: Vec<String>,
    pub hot_loop_nodes: Vec<String>,
    pub kinds: BTreeMap<String, usize>,
}

impl Default for Numerical {
    fn default() -> Self {
        let mut kinds = BTreeMap::from([
            ("iso_fortran_env::real32".into(), 32),
            ("iso_fortran_env::real64".into(), 64),
            ("iso_fortran_env::real128".into(), 128),
        ]);
        for module in [
            "iso_c_binding",
            "petscsys",
            "petscdmplex",
            "petscsnes",
            "petscvec",
        ] {
            for (name, width) in [
                ("c_float", 32),
                ("c_float_complex", 32),
                ("c_double", 64),
                ("c_double_complex", 64),
            ] {
                kinds.insert(format!("{module}::{name}"), width);
            }
        }
        Self {
            kernels: Vec::new(),
            hot_loops: Vec::new(),
            purity_warnings: Vec::new(),
            io_statements: list(&[
                "read_statement",
                "write_statement",
                "print_statement",
                "open_statement",
                "close_statement",
                "inquire_statement",
                "file_position_statement",
                "wait_statement",
            ]),
            hot_loop_nodes: list(&["do_loop", "while_statement"]),
            kinds,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Layer {
    pub paths: Vec<String>,
    pub may_use: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct External {
    pub layers: Vec<String>,
    pub paths: Vec<String>,
    pub intrinsic: bool,
    pub compiler: bool,
    pub symbols: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Internal {
    pub module: String,
    pub layers: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RuleSetting {
    pub severity: Option<Severity>,
    pub paths: Vec<String>,
    pub exclude_paths: Vec<String>,
    pub subjects: Vec<String>,
    pub exclude_subjects: Vec<String>,
    pub exception_categories: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RuleOverride {
    pub rules: Vec<String>,
    #[serde(default)]
    pub paths: Vec<String>,
    #[serde(default)]
    pub subjects: Vec<String>,
    pub severity: Severity,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Exception {
    pub rule: String,
    pub subject: String,
    pub category: String,
    pub reason: String,
    #[serde(default)]
    pub paths: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Inputs {
    pub extensions: Vec<String>,
    pub include: Vec<String>,
    pub exclude: Vec<String>,
    pub discovery_exclude: Vec<String>,
    pub follow_links: bool,
    pub allow_empty: bool,
}

impl Default for Inputs {
    fn default() -> Self {
        Self {
            extensions: list(&["f90", "f95", "f03", "f08", "f18", "fypp"]),
            include: Vec::new(),
            exclude: list(&["**/third_party/**"]),
            discovery_exclude: list(&[
                "**/.git/**",
                "**/.venv/**",
                "**/target/**",
                "**/.fortran-check/**",
                "**/build/**",
                "**/CMakeBuild*/**",
                "**/CMakeFiles/**",
            ]),
            follow_links: false,
            allow_empty: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Fortitude {
    pub enabled: bool,
    pub executable: PathBuf,
    pub version: String,
    pub config: Option<PathBuf>,
    pub select: Vec<String>,
    pub ignore: Vec<String>,
    pub preview: bool,
    pub allow_comments: bool,
    pub severity: Severity,
    pub batch_size: usize,
}

impl Default for Fortitude {
    fn default() -> Self {
        Self {
            enabled: true,
            executable: "fortitude".into(),
            version: "0.10.1".into(),
            config: None,
            select: list(&["E", "OB", "PORT021", "S001"]),
            ignore: Vec::new(),
            preview: false,
            allow_comments: false,
            severity: Severity::Error,
            batch_size: 128,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Compiler {
    pub gnu_flags: Vec<String>,
    pub ifx_flags: Vec<String>,
    pub ifx_windows_flags: Vec<String>,
    pub gnu_cpp: Vec<String>,
    pub ifx_cpp: Vec<String>,
    pub ifx_windows_cpp: Vec<String>,
}

impl Default for Compiler {
    fn default() -> Self {
        Self {
            gnu_flags: list(&[
                "-std={standard}",
                "-pedantic-errors",
                "-fimplicit-none",
                "-Wall",
                "-Wextra",
                "-Wimplicit-interface",
                "-Wimplicit-procedure",
                "-Werror",
            ]),
            ifx_flags: list(&[
                "-stand",
                "{intel_standard}",
                "-warn",
                "all,errors,stderrors,interfaces,externals",
                "-implicitnone",
            ]),
            ifx_windows_flags: list(&[
                "/stand:{intel_standard}",
                "/warn:all,errors,stderrors,interfaces,externals",
                "/implicitnone",
            ]),
            gnu_cpp: list(&["-cpp", "-E"]),
            ifx_cpp: list(&["-fpp", "-E"]),
            ifx_windows_cpp: list(&["/fpp", "/E"]),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Fypp {
    pub executable: PathBuf,
    pub args: Vec<String>,
    pub overrides: Vec<FyppOverride>,
}

impl Default for Fypp {
    fn default() -> Self {
        Self {
            executable: "fypp".into(),
            args: list(&[
                "-n",
                "--line-marker-format=std",
                "--line-length={line_length}",
            ]),
            overrides: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FyppOverride {
    pub paths: Vec<String>,
    pub args: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Runtime {
    pub jobs: usize,
    pub cache: bool,
    pub cache_directory: Option<PathBuf>,
}

impl Default for Runtime {
    fn default() -> Self {
        Self {
            jobs: 8,
            cache: true,
            cache_directory: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Reports {
    pub include_sarif: bool,
    pub include_metrics: bool,
    pub include_graph: bool,
}

impl Default for Reports {
    fn default() -> Self {
        Self {
            include_sarif: true,
            include_metrics: true,
            include_graph: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Config {
    pub version: u32,
    pub limits: Limits,
    pub metrics: Metrics,
    pub language: Language,
    pub numerical: Numerical,
    pub layers: BTreeMap<String, Layer>,
    pub externals: BTreeMap<String, External>,
    pub internals: Vec<Internal>,
    pub rules: BTreeMap<String, RuleSetting>,
    pub overrides: Vec<RuleOverride>,
    pub exceptions: Vec<Exception>,
    pub safety: SafetyPolicy,
    pub inputs: Inputs,
    pub fortitude: Fortitude,
    pub compiler: Compiler,
    pub fypp: Fypp,
    pub runtime: Runtime,
    pub reports: Reports,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            version: 1,
            limits: Limits::default(),
            metrics: Metrics::default(),
            language: Language::default(),
            numerical: Numerical::default(),
            layers: BTreeMap::new(),
            externals: BTreeMap::new(),
            internals: Vec::new(),
            rules: BTreeMap::new(),
            overrides: Vec::new(),
            exceptions: Vec::new(),
            safety: SafetyPolicy::default(),
            inputs: Inputs::default(),
            fortitude: Fortitude::default(),
            compiler: Compiler::default(),
            fypp: Fypp::default(),
            runtime: Runtime::default(),
            reports: Reports::default(),
        }
    }
}

impl Config {
    pub fn from_toml(text: &str) -> Result<Self> {
        let mut config: Self = toml::from_str(text).context("invalid policy TOML")?;
        config.validate()?;
        Ok(config)
    }

    pub fn load(profile: &str, overlay: Option<&Path>) -> Result<Self> {
        let contents = match profile {
            "generic" => include_str!("../profiles/generic.toml"),
            "ftcms" => include_str!("../profiles/ftcms.toml"),
            name => bail!("unknown built-in profile '{name}'"),
        };
        let mut value = toml::Value::try_from(Self::default())?;
        merge(
            &mut value,
            toml::from_str(include_str!("../profiles/fortran2018.toml"))?,
        );
        merge(&mut value, toml::from_str(contents)?);
        if let Some(path) = overlay {
            let raw = fs::read_to_string(path)
                .with_context(|| format!("cannot read policy {}", path.display()))?;
            merge(
                &mut value,
                toml::from_str(&raw)
                    .with_context(|| format!("invalid TOML in {}", path.display()))?,
            );
        }
        let mut config: Self = value
            .try_into()
            .context("unknown or invalid policy setting")?;
        config.validate()?;
        Ok(config)
    }

    pub fn filter(&self, mut diagnostic: Diagnostic, root: &Path) -> Option<Diagnostic> {
        let path = slash_path(
            diagnostic
                .location
                .path
                .strip_prefix(root)
                .unwrap_or(&diagnostic.location.path),
        );
        if self.exceptions.iter().any(|entry| {
            entry.rule == diagnostic.code
                && entry.subject.eq_ignore_ascii_case(&diagnostic.subject)
                && included(&entry.paths, &path)
        }) {
            return None;
        }
        if let Some(rule) = self.rules.get(&diagnostic.code) {
            if !included(&rule.paths, &path)
                || !included(&rule.subjects, &diagnostic.subject)
                || any_match(&rule.exclude_paths, &path)
                || any_match(&rule.exclude_subjects, &diagnostic.subject)
            {
                return None;
            }
            if let Some(severity) = rule.severity {
                diagnostic.severity = severity;
            }
        }
        for entry in &self.overrides {
            if any_match(&entry.rules, &diagnostic.code)
                && included(&entry.paths, &path)
                && included(&entry.subjects, &diagnostic.subject)
            {
                diagnostic.severity = entry.severity;
            }
        }
        (diagnostic.severity != Severity::Off).then_some(diagnostic)
    }

    fn validate(&mut self) -> Result<()> {
        ensure!(self.version == 1, "policy version must be 1");
        ensure!(self.runtime.jobs > 0, "runtime.jobs must be positive");
        ensure!(
            self.language.line_length > 0,
            "language.line_length must be positive"
        );
        ensure!(
            ["f95", "f2003", "f2008", "f2018", "f2023"].contains(&self.language.standard.as_str()),
            "unknown Fortran standard '{}', expected f95 through f2023",
            self.language.standard
        );
        ensure!(
            self.fortitude.batch_size > 0,
            "fortitude.batch_size must be positive"
        );
        ensure!(
            !self.fortitude.enabled || !self.fortitude.select.is_empty(),
            "fortitude.select must name one or more rules"
        );
        self.rules = normalize_keys(
            std::mem::take(&mut self.rules),
            str::to_ascii_uppercase,
            "rule",
        )?;
        self.externals = normalize_keys(
            std::mem::take(&mut self.externals),
            str::to_ascii_lowercase,
            "external module",
        )?;
        self.layers = normalize_keys(
            std::mem::take(&mut self.layers),
            str::to_ascii_lowercase,
            "layer",
        )?;
        self.limits.subjects = normalize_keys(
            std::mem::take(&mut self.limits.subjects),
            str::to_ascii_lowercase,
            "metric subject",
        )?;
        for (name, limits) in &self.limits.subjects {
            ensure!(qualified_subject(name), "invalid metric subject '{name}'");
            ensure!(
                !limits.reason.trim().is_empty(),
                "metric subject '{name}' requires a reason"
            );
            ensure!(
                limits.lcom4.is_some()
                    || limits.cbo.is_some()
                    || limits.wmc.is_some()
                    || limits.complexity.is_some(),
                "metric subject '{name}' requires a limit"
            );
        }
        self.numerical.kinds = normalize_keys(
            std::mem::take(&mut self.numerical.kinds),
            str::to_ascii_lowercase,
            "numeric kind",
        )?;
        for (name, layer) in &mut self.layers {
            ensure!(
                !layer.paths.is_empty(),
                "layers.{name} requires one or more paths"
            );
            validate_globs(&layer.paths)?;
            layer
                .may_use
                .iter_mut()
                .for_each(|item| item.make_ascii_lowercase());
        }
        for (name, external) in &mut self.externals {
            ensure!(
                !external.compiler || external.symbols.is_empty(),
                "externals.{name} cannot combine a symbol catalog with compiler resolution"
            );
            external
                .layers
                .iter_mut()
                .for_each(|item| item.make_ascii_lowercase());
            external
                .symbols
                .iter_mut()
                .for_each(|item| item.make_ascii_lowercase());
            validate_globs(&external.paths)?;
        }
        for internal in &mut self.internals {
            internal
                .layers
                .iter_mut()
                .for_each(|item| item.make_ascii_lowercase());
            validate_globs(std::slice::from_ref(&internal.module))?;
        }
        for name in self
            .layers
            .values()
            .flat_map(|layer| &layer.may_use)
            .chain(self.externals.values().flat_map(|module| &module.layers))
            .chain(self.internals.iter().flat_map(|rule| &rule.layers))
        {
            ensure!(
                name == "*" || self.layers.contains_key(name),
                "unknown architecture layer '{name}'"
            );
        }
        for (name, rule) in &self.rules {
            ensure!(valid_rule(name), "unknown rule ID '{name}'");
            for patterns in [
                &rule.paths,
                &rule.exclude_paths,
                &rule.subjects,
                &rule.exclude_subjects,
            ] {
                validate_globs(patterns)?;
            }
            if let Some(categories) = &rule.exception_categories {
                ensure!(
                    name != "STD002"
                        || categories
                            .iter()
                            .all(|category| category == "significant_benefit"),
                    "INTENT(OUT) exceptions require significant_benefit"
                );
                ensure!(
                    categories
                        .iter()
                        .all(|category| category == "*" || !category.trim().is_empty()),
                    "empty exception category for {name}"
                );
            }
        }
        for entry in &self.overrides {
            ensure!(
                !entry.rules.is_empty(),
                "an override requires one or more rule selectors"
            );
            for selector in &entry.rules {
                ensure!(
                    selector.contains(['*', '?', '[', '{'])
                        || valid_rule(&selector.to_ascii_uppercase()),
                    "unknown override rule '{selector}'"
                );
            }
            for patterns in [&entry.rules, &entry.paths, &entry.subjects] {
                validate_globs(patterns)?;
            }
        }
        let mut seen = BTreeSet::new();
        for entry in &mut self.exceptions {
            entry.rule.make_ascii_uppercase();
            entry.subject.make_ascii_lowercase();
            entry.category.make_ascii_lowercase();
            ensure!(
                valid_rule(&entry.rule),
                "unknown exception rule '{}', cannot apply exception",
                entry.rule
            );
            ensure!(
                qualified_subject(&entry.subject),
                "exceptions require an exact qualified subject"
            );
            ensure!(
                !entry.reason.trim().is_empty() && !entry.category.trim().is_empty(),
                "an exception requires a category and reason"
            );
            let configured = self
                .rules
                .get(&entry.rule)
                .and_then(|rule| rule.exception_categories.as_ref());
            let permitted = configured.map_or_else(
                || default_categories(&entry.rule).contains(&entry.category.as_str()),
                |categories| {
                    categories.iter().any(|category| {
                        category == "*" || category.eq_ignore_ascii_case(&entry.category)
                    })
                },
            );
            ensure!(
                permitted,
                "exception category '{}' is unavailable for {}",
                entry.category,
                entry.rule
            );
            validate_globs(&entry.paths)?;
            ensure!(
                seen.insert((&entry.rule, &entry.subject, &entry.paths)),
                "duplicate exception for {} on {}",
                entry.rule,
                entry.subject
            );
        }
        for patterns in [
            &self.inputs.include,
            &self.inputs.exclude,
            &self.inputs.discovery_exclude,
            &self.numerical.kernels,
            &self.numerical.hot_loops,
            &self.numerical.purity_warnings,
        ] {
            validate_globs(patterns)?;
        }
        for rule in &self.fypp.overrides {
            ensure!(
                !rule.paths.is_empty(),
                "a Fypp override requires one or more paths"
            );
            validate_globs(&rule.paths)?;
        }
        self.safety.validate()?;
        Ok(())
    }
}

fn normalize_keys<T>(
    input: BTreeMap<String, T>,
    normalize: impl Fn(&str) -> String,
    label: &str,
) -> Result<BTreeMap<String, T>> {
    let mut result = BTreeMap::new();
    for (name, value) in input {
        let key = normalize(&name);
        ensure!(!result.contains_key(&key), "duplicate {label} '{name}'");
        result.insert(key, value);
    }
    Ok(result)
}

fn merge(base: &mut toml::Value, overlay: toml::Value) {
    match (base, overlay) {
        (toml::Value::Table(base_table), toml::Value::Table(overlay_table)) => {
            for (name, value) in overlay_table {
                if let Some(existing) = base_table.get_mut(&name) {
                    merge(existing, value);
                } else {
                    base_table.insert(name, value);
                }
            }
        }
        // Scalar and array settings are replaced together.
        (base, overlay) => *base = overlay,
    }
}

pub fn list(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| (*item).into()).collect()
}

pub fn matches(pattern: &str, value: &str) -> bool {
    if pattern == "*" {
        return true;
    }
    if !pattern.contains(['*', '?', '[', '{']) {
        return pattern.eq_ignore_ascii_case(value);
    }
    thread_local! {
        static MATCHERS: RefCell<HashMap<String, Option<GlobMatcher>>> = RefCell::new(HashMap::new());
    }
    MATCHERS.with_borrow_mut(|matchers| {
        if let Some(matcher) = matchers.get(pattern) {
            return matcher
                .as_ref()
                .is_some_and(|matcher| matcher.is_match(value));
        }
        let matcher = glob(pattern)
            .build()
            .ok()
            .map(|glob| glob.compile_matcher());
        let result = matcher
            .as_ref()
            .is_some_and(|matcher| matcher.is_match(value));
        matchers.insert(pattern.to_owned(), matcher);
        result
    })
}

pub fn any_match(patterns: &[String], value: &str) -> bool {
    patterns.iter().any(|pattern| matches(pattern, value))
}

pub fn included(patterns: &[String], value: &str) -> bool {
    patterns.is_empty() || any_match(patterns, value)
}

fn validate_globs(patterns: &[String]) -> Result<()> {
    for pattern in patterns {
        ensure!(!pattern.is_empty(), "empty glob selector");
        glob(pattern)
            .build()
            .with_context(|| format!("invalid glob '{pattern}'"))?;
    }
    Ok(())
}

fn glob(pattern: &str) -> GlobBuilder<'_> {
    let mut builder = GlobBuilder::new(pattern);
    builder.case_insensitive(true).literal_separator(true);
    builder
}

fn qualified_subject(subject: &str) -> bool {
    static SUBJECT: OnceLock<Regex> = OnceLock::new();
    SUBJECT.get_or_init(|| {
        let entity = r"(?:[a-z][a-z0-9_]*(?::[a-z][a-z0-9_]*)?|(?:assignment|operator|read|write)\([.a-z0-9_+\-*/<>=]+\))";
        Regex::new(&format!(r"\A{entity}(?:::{entity})*\z")).expect("valid qualified Fortran name pattern")
    }).is_match(subject)
}

fn default_categories(rule: &str) -> &[&str] {
    match rule {
        "STD002" => &["significant_benefit"],
        "OOP002" => &["interop"],
        "OOP005" | "OOP007" => &["constructor", "factory", "callback", "interop"],
        "OOP008" => &["polymorphism", "callback", "interop"],
        "OOP009" => &["callback", "interop", "method_without_object"],
        "OOP010" => &["borrowed_pointer", "interop"],
        "COH001" => &["cohesion"],
        _ => &[],
    }
}

fn valid_rule(id: &str) -> bool {
    RULES.iter().any(|rule| rule.id == id)
        || ["C", "E", "OB", "PORT", "S", "MOD", "IO", "P", "M", "R"]
            .iter()
            .any(|prefix| {
                id.strip_prefix(prefix).is_some_and(|suffix| {
                    suffix.len() >= 3 && suffix.bytes().all(|byte| byte.is_ascii_digit())
                })
            })
}
