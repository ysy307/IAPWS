use std::path::Path;

use fortran_check::config::{Config, External, Internal, Layer};
use fortran_check::design;
use fortran_check::semantic::Project;
use fortran_check::source::{Source, Unit};

fn analyze(inputs: &[(&str, &str)], config: &Config) -> design::DesignReport {
    let units: Vec<_> = inputs
        .iter()
        .map(|(path, text)| Unit::parse(Source::new(*path, *text)).unwrap())
        .collect();
    let project = Project::new(&units, config);
    design::analyze(&project, config, Path::new("."))
}

fn codes(report: &design::DesignReport) -> Vec<&str> {
    report
        .diagnostics
        .iter()
        .map(|item| item.code.as_str())
        .collect()
}

#[test]
fn inheritance_counts_edges_children_and_unresolved_depths() {
    let mut config = Config::default();
    config.limits.dit = 1;
    let report = analyze(
        &[(
            "types.f90",
            "module inheritance
type :: base
end type
type, extends(base) :: middle
end type
type, extends(base) :: sibling
end type
type, extends(middle) :: leaf
end type
type, extends(missing) :: orphan
end type
type, extends(cycle_b) :: cycle_a
end type
type, extends(cycle_a) :: cycle_b
end type
end module",
        )],
        &config,
    );
    let types = &report.metrics["types"];
    assert_eq!(types["inheritance::base"]["inheritance"]["dit"], 0);
    assert_eq!(types["inheritance::base"]["inheritance"]["noc"], 2);
    assert_eq!(types["inheritance::leaf"]["inheritance"]["dit"], 2);
    assert!(types["inheritance::orphan"]["inheritance"]["dit"].is_null());
    assert!(types["inheritance::cycle_a"]["inheritance"]["dit"].is_null());
    for code in ["DIT001", "DIT002", "DIT003"] {
        assert!(codes(&report).contains(&code), "missing {code}");
    }
}

#[test]
fn cohesion_and_coupling_use_resolved_methods_and_passed_object_components() {
    let mut config = Config::default();
    config.limits.wmc = 3;
    config.limits.cbo = 0;
    config.limits.complexity = 1;
    let report = analyze(
        &[(
            "model.f90",
            "module model
type :: dependency
end type
type :: counter
   integer :: left, right
   type(dependency) :: state
contains
   procedure :: a, b, read_right, relay
end type
contains
subroutine a(self)
   class(counter) :: self
   if (self%left > 0) self%left = 0
end subroutine
subroutine b(self)
   class(counter) :: self
   self%left = self%left + 1
end subroutine
integer function read_right(self) result(value)
   class(counter) :: self
   value = self%right
end function
integer function relay(self) result(value)
   class(counter) :: self
   value = self%read_right()
end function
end module",
        )],
        &config,
    );
    let metrics = &report.metrics["types"]["model::counter"];
    assert_eq!(metrics["cohesion"]["lcom4"], 2);
    assert_eq!(metrics["coupling"]["cbo"], 1);
    assert_eq!(metrics["coupling"]["wmc"], 5);
    assert_eq!(metrics["coupling"]["rfc"], 4);
    for code in ["COH001", "CBO001", "WMC001", "CC001"] {
        assert!(codes(&report).contains(&code), "missing {code}");
    }
}

#[test]
fn architecture_resolves_layers_reports_cycles_and_uses_distinct_edges() {
    let mut config = Config::default();
    config.layers.clear();
    config.layers.insert(
        "core".into(),
        Layer {
            paths: vec!["core/**".into()],
            may_use: vec!["core".into()],
        },
    );
    config.layers.insert(
        "app".into(),
        Layer {
            paths: vec!["app/**".into()],
            may_use: vec!["core".into(), "app".into()],
        },
    );
    config.externals.insert(
        "vendor".into(),
        External {
            layers: vec!["app".into()],
            paths: vec!["app/**".into()],
            compiler: true,
            ..Default::default()
        },
    );
    config.internals.push(Internal {
        module: "inner".into(),
        layers: vec!["core".into()],
    });
    let inputs = [
        (
            "core/inner.f90",
            "module inner\nuse runner\nuse vendor\nend module",
        ),
        (
            "app/main.f90",
            "module runner\nuse inner\nuse inner\nuse missing\nend module",
        ),
        ("app/second.f90", "module Runner\nend module"),
        ("other.f90", "module orphan\nend module"),
    ];
    let report = analyze(&inputs, &config);
    for code in [
        "ARCH001", "ARCH002", "ARCH003", "ARCH004", "ARCH005", "ARCH006", "ARCH007",
    ] {
        assert!(codes(&report).contains(&code), "missing {code}");
    }
    assert_eq!(report.metrics["modules"]["inner"]["fan_in"], 1);
    assert_eq!(report.metrics["modules"]["inner"]["fan_out"], 2);
    assert_eq!(report.metrics["modules"]["runner"]["fan_out"], 1);
    assert_eq!(report.metrics["modules"]["orphan"]["instability"], 0.0);
    assert!(report.graph.contains("\"inner\" -> \"external::vendor\";"));
    assert!(report.graph.contains("\"runner\" -> \"inner\";"));
    config.metrics.external_fan_out = false;
    let internal_only = analyze(&inputs, &config);
    assert_eq!(internal_only.metrics["modules"]["inner"]["fan_out"], 1);
    assert_eq!(internal_only.graph, report.graph);
}

#[test]
fn architecture_recognizes_intrinsic_shadowing_and_external_reexports() {
    let mut config = Config::load("generic", None).unwrap();
    config.externals.insert(
        "vendor".into(),
        External {
            paths: vec!["bridge/**".into()],
            compiler: true,
            ..Default::default()
        },
    );
    let report = analyze(
        &[
            (
                "bridge/bridge.f90",
                "module bridge\nuse vendor, only: solve => vendor_solve\nprivate\npublic :: solve\nend module",
            ),
            (
                "src/client.f90",
                "module client\nuse bridge, only: solve\nend module",
            ),
            ("src/private_env.f90", "module iso_fortran_env\nend module"),
            (
                "src/mixed.f90",
                "module mixed\nuse, intrinsic :: iso_fortran_env, only: real64\nuse, non_intrinsic :: iso_fortran_env\nend module",
            ),
        ],
        &config,
    );
    assert!(
        report
            .graph
            .contains("\"mixed\" -> \"external::iso_fortran_env\";")
    );
    assert!(report.graph.contains("\"mixed\" -> \"iso_fortran_env\";"));
    assert!(report.graph.contains("\"client\" -> \"external::vendor\";"));
    assert!(
        report
            .diagnostics
            .iter()
            .any(|item| item.code == "ARCH006" && item.subject == "client::vendor")
    );
}

#[test]
fn cohesion_exemptions_and_complexity_nodes_are_configurable() {
    let source = r#"module model
type :: stateless
contains
    procedure :: first
    procedure :: second
end type
type, abstract :: prototype
    integer :: left, right
contains
    procedure :: reset_left
    procedure :: reset_right
end type
contains
subroutine first(self)
    class(stateless) :: self
    if (.true.) return
end subroutine
subroutine second(self)
    class(stateless) :: self
end subroutine
subroutine reset_left(self)
    class(prototype) :: self
    self%left = 0
end subroutine
subroutine reset_right(self)
    class(prototype) :: self
    self%right = 0
end subroutine
end module
"#;
    let mut config = Config::default();
    let exempt = analyze(&[("model.f90", source)], &config);
    assert!(exempt.diagnostics.is_empty(), "{:?}", exempt.diagnostics);
    assert_eq!(
        exempt.metrics["types"]["model::stateless"]["cohesion"]["lcom4"],
        2
    );
    assert_eq!(
        exempt.metrics["types"]["model::stateless"]["coupling"]["wmc"],
        3
    );
    config.metrics.exclude_abstract = false;
    config.metrics.exclude_stateless = false;
    config.metrics.complexity_nodes.clear();
    let strict = analyze(&[("model.f90", source)], &config);
    assert_eq!(codes(&strict), vec!["COH001", "COH001"]);
    assert_eq!(
        strict.metrics["types"]["model::stateless"]["coupling"]["wmc"],
        2
    );
}

#[test]
fn qualified_metric_limits_preserve_measurements_and_flag_growth() {
    let source = r#"
module model
type :: dependency
end type
type :: counter
    integer :: left, right
    type(dependency) :: state
contains
    procedure :: update, reset
end type
contains
subroutine update(self)
    class(counter) :: self
    if (self%left > 0) self%left = 0
end subroutine
subroutine reset(self)
    class(counter) :: self
    self%right = 0
end subroutine
subroutine grows(value)
    integer :: value
    if (value > 1) value = 1
    if (value < 0) value = 0
end subroutine
subroutine unrelated(value)
    integer :: value
    if (value > 0) value = 0
end subroutine
end module
"#;
    let config = Config::from_toml(
        r#"
[limits]
cbo = 0
wmc = 2
complexity = 1

[limits.subjects."Model::Counter"]
lcom4 = 2
cbo = 1
wmc = 3
reason = "The published bindings work on distinct fields."

[limits.subjects."model::counter::update"]
complexity = 2
reason = "The update uses a range guard."

[limits.subjects."model::grows"]
complexity = 2
reason = "The first boundary check needs one branch."
"#,
    )
    .unwrap();
    let report = analyze(&[("model.f90", source)], &config);
    assert_eq!(
        report.metrics,
        analyze(&[("model.f90", source)], &Config::default()).metrics
    );
    let violations: Vec<_> = report
        .diagnostics
        .iter()
        .map(|item| (item.code.as_str(), item.subject.as_str()))
        .collect();
    assert_eq!(
        violations,
        [("CC001", "model::grows"), ("CC001", "model::unrelated")]
    );
}
