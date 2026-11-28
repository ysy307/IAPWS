//! Index lexical declarations and resolve references across translation units.

use crate::config::Config;
use crate::diagnostic::{Diagnostic, Location};
use crate::source::{Unit, child, children, fields};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use tree_sitter::Node;

pub const PROGRAM_UNITS: &[&str] = &[
    "module",
    "submodule",
    "program",
    "subroutine",
    "function",
    "module_procedure",
];
const NAME_NODES: &[&str] = &[
    "name",
    "identifier",
    "method_name",
    "type_name",
    "operator",
    "assignment",
    "defined_io_procedure",
];
const CALLABLE: &[&str] = &["procedure", "generic", "binding"];
const TYPES: &[&str] = &["type"];
const PROCEDURES: &[&str] = &["procedure"];
const INTERFACES: &[&str] = &["procedure", "abstract_interface", "generic"];
const SKIP: &[&str] = &[
    "use_statement",
    "import_statement",
    "implicit_statement",
    "internal_procedures",
    "interface",
    "comment",
    "statement_label",
    "language_binding",
    "construct_name",
    "construct_name_label",
    "contains_statement",
    "procedure_attribute",
    "type_qualifier",
    "procedure_qualifier",
];
const INTRINSICS: &str = "
abs achar acos acosh adjustl adjustr aimag aint all allocated anint any asin asinh associated atan atan2 atanh
atomic_add atomic_and atomic_cas atomic_define atomic_fetch_add atomic_fetch_and atomic_fetch_or atomic_fetch_xor
atomic_or atomic_ref atomic_xor bessel_j0 bessel_j1 bessel_jn bessel_y0 bessel_y1 bessel_yn bge bgt bit_size ble blt
btest ceiling char cmplx co_broadcast co_max co_min co_reduce co_sum command_argument_count conjg cos cosh count
cpu_time cshift date_and_time dble digits dim dot_product dprod dshiftl dshiftr eoshift epsilon erf erfc erfc_scaled
event_query execute_command_line exp exponent extends_type_of failed_images findloc floor fraction gamma
get_command get_command_argument get_environment_variable get_team huge hypot iachar iall iand iany ibclr ibits
ibset ichar ieor image_index image_status index int ior iparity is_contiguous is_iostat_end is_iostat_eor ishft
ishftc kind lbound lcobound leadz len len_trim lge lgt lle llt log log10 log_gamma logical maskl maskr matmul max
maxexponent maxloc maxval merge merge_bits min minexponent minloc minval mod modulo move_alloc mvbits nearest
new_line nint norm2 not null num_images out_of_range pack parity popcnt poppar precision present product radix
random_init random_number random_seed range rank real reduce repeat reshape rrspacing same_type_as scale scan
selected_char_kind selected_int_kind selected_real_kind set_exponent shape shifta shiftl shiftr sign sin sinh size
spacing spread sqrt stopped_images storage_size sum system_clock tan tanh team_number this_image tiny trailz
transfer transpose trim ubound ucobound unpack verify
";

#[derive(Debug, Clone)]
pub struct Use {
    pub module: String,
    pub location: Location,
    pub intrinsic: bool,
    pub non_intrinsic: bool,
    pub only: bool,
    pub imports: BTreeMap<String, String>,
    pub imported_at: BTreeMap<String, Location>,
}

impl Use {
    fn remote<'a>(&'a self, local: &'a str) -> Option<&'a str> {
        if let Some(remote) = self.imports.get(local) {
            Some(remote)
        } else if self.only || self.imports.values().any(|remote| remote == local) {
            None
        } else {
            Some(local)
        }
    }
}

#[derive(Debug)]
pub struct Scope<'a> {
    pub kind: &'static str,
    pub name: String,
    pub qualified: String,
    pub unit: usize,
    pub node: Node<'a>,
    pub header: Node<'a>,
    pub parent: Option<usize>,
    pub attrs: BTreeSet<String>,
    pub interface_body: bool,
    pub parameters: Vec<String>,
    pub result: Option<String>,
    pub symbols: BTreeMap<String, usize>,
    pub generics: BTreeMap<String, usize>,
    pub uses: Vec<Use>,
    pub default_private: bool,
    pub declaration: Option<usize>,
    pub ancestor: Option<String>,
    pub parent_name: Option<String>,
    pub explicit_access: BTreeMap<String, bool>,
    pub explicit_access_nodes: BTreeMap<String, Node<'a>>,
    imports: Option<BTreeSet<String>>,
    symbol: Option<usize>,
}

impl Scope<'_> {
    fn public(&self, name: &str, attrs: &BTreeSet<String>) -> bool {
        if attrs.contains("private") {
            false
        } else if attrs.contains("public") {
            true
        } else {
            self.explicit_access
                .get(name)
                .copied()
                .unwrap_or(!self.default_private)
        }
    }
}

#[derive(Debug)]
pub struct Symbol<'a> {
    pub name: String,
    pub qualified: String,
    pub kind: &'static str,
    pub scope: usize,
    pub node: Node<'a>,
    pub attrs: BTreeSet<String>,
    pub is_public: bool,
    pub type_name: Option<String>,
    pub intrinsic: Option<String>,
    pub target: Option<String>,
    pub base_name: Option<String>,
    pub members: BTreeMap<String, usize>,
    pub specifics: Vec<String>,
    pub owner: Option<usize>,
    pub procedure: Option<usize>,
    pub external: Option<String>,
    pub interface_name: Option<String>,
    rank: Option<usize>,
}

#[derive(Debug, Clone, Copy)]
pub struct Reference<'a> {
    pub scope: usize,
    pub node: Node<'a>,
    pub symbol: usize,
    pub receiver: Option<usize>,
    pub root_receiver: Option<usize>,
    pub call: bool,
}

#[derive(Debug)]
struct ExternalNames {
    intrinsic: bool,
    compiler: bool,
    symbols: BTreeSet<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum Candidate {
    Local(usize),
    External(String, String),
}

type Aliases = BTreeMap<String, Option<usize>>;

struct Argument {
    keyword: Option<String>,
    category: Option<String>,
    derived: Option<usize>,
    rank: Option<usize>,
}

pub struct Project<'a> {
    pub units: &'a [Unit],
    pub scopes: Vec<Scope<'a>>,
    pub symbols: Vec<Symbol<'a>>,
    pub modules: BTreeMap<String, usize>,
    pub types: Vec<usize>,
    pub references: Vec<Reference<'a>>,
    pub diagnostics: Vec<Diagnostic>,
    pub external_imports: Vec<(usize, usize, Location)>,
    externals: BTreeMap<String, ExternalNames>,
    external_symbols: BTreeMap<(String, String), usize>,
    external_bodies: BTreeMap<String, usize>,
    intrinsics: BTreeSet<&'static str>,
    submodules: BTreeMap<String, usize>,
    hosts: BTreeMap<usize, usize>,
    implementations: BTreeMap<usize, usize>,
    scope_by_node: HashMap<(usize, usize), usize>,
    expression_symbols: HashMap<(usize, usize), usize>,
    receiver_roots: HashMap<(usize, usize), usize>,
}

impl<'a> Project<'a> {
    pub fn new(units: &'a [Unit], config: &Config) -> Self {
        let externals = config
            .externals
            .iter()
            .map(|(name, external)| {
                (
                    name.to_ascii_lowercase(),
                    ExternalNames {
                        intrinsic: external.intrinsic,
                        compiler: external.compiler,
                        symbols: external
                            .symbols
                            .iter()
                            .map(|name| name.to_ascii_lowercase())
                            .collect(),
                    },
                )
            })
            .collect();
        let mut result = Self {
            units,
            scopes: Vec::new(),
            symbols: Vec::new(),
            modules: BTreeMap::new(),
            types: Vec::new(),
            references: Vec::new(),
            diagnostics: Vec::new(),
            external_imports: Vec::new(),
            externals,
            external_symbols: BTreeMap::new(),
            external_bodies: BTreeMap::new(),
            intrinsics: INTRINSICS.split_whitespace().collect(),
            submodules: BTreeMap::new(),
            hosts: BTreeMap::new(),
            implementations: BTreeMap::new(),
            scope_by_node: HashMap::new(),
            expression_symbols: HashMap::new(),
            receiver_roots: HashMap::new(),
        };
        for (index, unit) in units.iter().enumerate() {
            for node in children(unit.tree.root_node()) {
                if PROGRAM_UNITS.contains(&node.kind()) && !node.has_error() {
                    result.add_scope(index, node, None, false, false);
                }
            }
        }
        result.link_implementations();
        result.check_declarations();
        for scope in 0..result.scopes.len() {
            if result.scopes[scope].kind != "block" {
                result.scan_body(scope, &Aliases::new());
            }
        }
        result
    }

    pub fn unit(&self, scope: usize) -> &'a Unit {
        &self.units[self.scopes[scope].unit]
    }

    pub fn host(&self, scope: usize) -> Option<usize> {
        self.scopes[scope]
            .parent
            .or_else(|| self.hosts.get(&scope).copied())
    }

    pub fn scope_of(&self, unit: usize, node: Node<'_>) -> Option<usize> {
        let mut next = Some(node);
        while let Some(current) = next {
            if let Some(scope) = self.scope_by_node.get(&(unit, current.id())) {
                return Some(*scope);
            }
            next = current.parent();
        }
        None
    }

    fn owner_scope(&self, mut scope: usize) -> usize {
        let mut seen = HashSet::new();
        while seen.insert(scope) {
            let Some(host) = self.host(scope) else { break };
            scope = host;
        }
        scope
    }

    fn issue(&mut self, code: &str, scope: usize, node: Node<'_>, name: &str, message: String) {
        let subject = if name.is_empty() {
            self.scopes[scope].qualified.clone()
        } else {
            format!("{}::{name}", self.scopes[scope].qualified)
        };
        self.diagnostics.push(
            Diagnostic::new(code, message, self.unit(scope).location(node)).with_subject(subject),
        );
    }

    fn symbol(&self, scope: usize, name: String, kind: &'static str, node: Node<'a>) -> Symbol<'a> {
        Symbol {
            qualified: format!("{}::{name}", self.scopes[scope].qualified),
            name,
            kind,
            scope,
            node,
            attrs: BTreeSet::new(),
            is_public: false,
            type_name: None,
            intrinsic: None,
            target: None,
            base_name: None,
            members: BTreeMap::new(),
            specifics: Vec::new(),
            owner: None,
            procedure: None,
            external: None,
            interface_name: None,
            rank: None,
        }
    }

    fn push_symbol(&mut self, symbol: Symbol<'a>) -> usize {
        let id = self.symbols.len();
        self.symbols.push(symbol);
        id
    }

    fn register(&mut self, mut symbol: Symbol<'a>, owner: Option<usize>) -> usize {
        let scope = symbol.scope;
        let node = symbol.node;
        let name = symbol.name.clone();
        if let Some(owner) = owner {
            symbol.qualified = format!("{}::{name}", self.symbols[owner].qualified);
            symbol.owner = Some(owner);
        }
        let duplicate = match owner {
            Some(owner) => self.symbols[owner].members.get(&name).copied(),
            None => self.scopes[scope].symbols.get(&name).copied(),
        };
        if let Some(existing) = duplicate {
            let subject = symbol.qualified;
            self.diagnostics.push(
                Diagnostic::new(
                    "SEM001",
                    format!("Duplicate declaration of '{name}'"),
                    self.unit(scope).location(node),
                )
                .with_subject(subject),
            );
            return existing;
        }
        let id = self.push_symbol(symbol);
        self.expression_symbols
            .insert((self.scopes[scope].unit, node.id()), id);
        match owner {
            Some(owner) => {
                self.symbols[owner].members.insert(name, id);
            }
            None => {
                self.scopes[scope].symbols.insert(name, id);
            }
        }
        id
    }

    fn add_scope(
        &mut self,
        unit_index: usize,
        node: Node<'a>,
        parent: Option<usize>,
        interface: bool,
        abstract_: bool,
    ) -> usize {
        let unit = &self.units[unit_index];
        let header = child(node, &format!("{}_statement", node.kind())).unwrap_or(node);
        let name_node = header
            .child_by_field_name("name")
            .or_else(|| child(header, "name"));
        let mut name = name_node
            .map(|part| unit.name(part))
            .unwrap_or_else(|| "<main>".into());
        let kind = if node.kind() == "block_construct" {
            name = format!("<block@{}>", node.start_position().row + 1);
            "block"
        } else {
            node.kind()
        };
        let ancestor = header
            .child_by_field_name("ancestor")
            .map(|part| unit.name(part));
        let parent_name = header
            .child_by_field_name("parent")
            .map(|part| unit.name(part));
        if kind == "submodule"
            && let Some(ancestor) = &ancestor
        {
            name = format!("{ancestor}:{name}");
        }
        let qualified = parent.map_or_else(
            || name.clone(),
            |id| format!("{}::{name}", self.scopes[id].qualified),
        );
        let attrs = attributes(unit, header);
        let parameters = header
            .child_by_field_name("parameters")
            .map(|part| {
                children(part)
                    .into_iter()
                    .filter(|item| item.kind() == "identifier")
                    .map(|item| unit.name(item))
                    .collect()
            })
            .unwrap_or_default();
        let result = (kind == "function").then(|| {
            child(header, "function_result")
                .and_then(|item| item.named_child(0))
                .map(|item| unit.name(item))
                .unwrap_or_else(|| name.clone())
        });
        let imports = (interface && !attrs.contains("module")).then(BTreeSet::new);
        let scope = self.scopes.len();
        self.scopes.push(Scope {
            kind,
            name: name.clone(),
            qualified,
            unit: unit_index,
            node,
            header,
            parent,
            attrs: attrs.clone(),
            interface_body: interface,
            parameters,
            result,
            symbols: BTreeMap::new(),
            generics: BTreeMap::new(),
            uses: Vec::new(),
            default_private: false,
            declaration: None,
            ancestor,
            parent_name,
            explicit_access: BTreeMap::new(),
            explicit_access_nodes: BTreeMap::new(),
            imports,
            symbol: None,
        });
        self.scope_by_node.insert((unit_index, node.id()), scope);
        match kind {
            "module" => {
                self.modules.entry(name).or_insert(scope);
            }
            "submodule" => {
                self.submodules.entry(name).or_insert(scope);
            }
            "subroutine" | "function" | "module_procedure" => {
                if let (Some(parent), Some(name_node)) = (parent, name_node) {
                    let mut symbol = self.symbol(
                        parent,
                        name.clone(),
                        if abstract_ {
                            "abstract_interface"
                        } else {
                            "procedure"
                        },
                        name_node,
                    );
                    symbol.is_public = self.scopes[parent].public(&name, &attrs);
                    symbol.attrs = attrs;
                    symbol.procedure = Some(scope);
                    let paired = !interface
                        && (kind == "module_procedure" || symbol.attrs.contains("module"))
                        && self
                            .local(parent, &name, PROCEDURES)
                            .and_then(|id| self.symbols[id].procedure)
                            .is_some_and(|id| {
                                self.scopes[id].interface_body
                                    && self.scopes[id].attrs.contains("module")
                            });
                    self.scopes[scope].symbol = Some(if paired {
                        self.push_symbol(symbol)
                    } else {
                        self.register(symbol, None)
                    });
                } else if parent.is_none() {
                    match self.external_bodies.entry(name) {
                        std::collections::btree_map::Entry::Vacant(entry) => {
                            entry.insert(scope);
                        }
                        std::collections::btree_map::Entry::Occupied(entry) => {
                            let message = format!("Duplicate external procedure '{}'", entry.key());
                            self.issue("SEM001", scope, header, "", message);
                        }
                    }
                }
            }
            _ => {}
        }
        self.specification(scope);
        if let Some(result) = self.scopes[scope].result.clone() {
            if !self.scopes[scope].symbols.contains_key(&result)
                && let Some(dtype) = children(header)
                    .into_iter()
                    .find(|item| matches!(item.kind(), "intrinsic_type" | "derived_type"))
            {
                let mut symbol = self.symbol(
                    scope,
                    result.clone(),
                    "variable",
                    name_node.unwrap_or(header),
                );
                symbol.type_name = declared_type(unit, Some(dtype));
                symbol.intrinsic = intrinsic_type(Some(dtype)).map(str::to_owned);
                symbol.rank = Some(0);
                self.register(symbol, None);
            }
            if let (Some(procedure), Some(entity)) = (
                self.scopes[scope].symbol,
                self.scopes[scope].symbols.get(&result).copied(),
            ) {
                self.symbols[procedure].type_name = self.symbols[entity].type_name.clone();
                self.symbols[procedure].intrinsic = self.symbols[entity].intrinsic.clone();
            }
        }
        scope
    }

    fn specification(&mut self, scope: usize) {
        let unit = self.unit(scope);
        let nodes = children(self.scopes[scope].node);
        for node in &nodes {
            if matches!(node.kind(), "private_statement" | "public_statement") {
                let public = node.kind() == "public_statement";
                let names: Vec<_> = children(*node)
                    .into_iter()
                    .filter(|part| NAME_NODES.contains(&part.kind()))
                    .collect();
                if names.is_empty() {
                    self.scopes[scope].default_private = !public;
                }
                for part in names {
                    let name = unit.name(part);
                    self.scopes[scope]
                        .explicit_access
                        .insert(name.clone(), public);
                    self.scopes[scope].explicit_access_nodes.insert(name, part);
                }
            }
        }
        for &node in &nodes {
            match node.kind() {
                "use_statement" => {
                    if let Some(use_) = use_statement(unit, node) {
                        self.scopes[scope].uses.push(use_);
                    }
                }
                "import_statement" => {
                    let names: BTreeSet<_> = children(node)
                        .into_iter()
                        .filter(|item| item.kind() == "identifier")
                        .map(|item| unit.name(item))
                        .collect();
                    let mode = spelling(unit, node);
                    if mode.starts_with("import,none") {
                        self.scopes[scope].imports = Some(BTreeSet::new());
                    } else if mode.starts_with("import,only") {
                        self.scopes[scope]
                            .imports
                            .get_or_insert_with(BTreeSet::new)
                            .extend(names);
                    } else if names.is_empty() {
                        self.scopes[scope].imports = None;
                    } else if let Some(existing) = &mut self.scopes[scope].imports {
                        existing.extend(names);
                    }
                }
                "variable_declaration" => self.variables(scope, node, None, false),
                "derived_type_definition" => self.derived_type(scope, node),
                "interface" => self.interface(scope, node),
                "internal_procedures" => {
                    for procedure in children(node)
                        .into_iter()
                        .filter(|item| PROGRAM_UNITS.contains(&item.kind()))
                    {
                        self.add_scope(
                            self.scopes[scope].unit,
                            procedure,
                            Some(scope),
                            false,
                            false,
                        );
                    }
                }
                _ => self.collect_blocks(scope, node),
            }
        }
        for node in nodes {
            match node.kind() {
                "variable_modification" => self.modify_variables(scope, node),
                "parameter_statement" => {
                    for part in children(node)
                        .into_iter()
                        .filter(|item| item.kind() == "parameter_assignment")
                    {
                        if let Some(name) = part.named_child(0).map(|item| unit.name(item))
                            && let Some(&id) = self.scopes[scope].symbols.get(&name)
                        {
                            self.symbols[id].kind = "constant";
                            self.symbols[id].attrs.insert("parameter".into());
                        }
                    }
                }
                _ => {}
            }
        }
    }

    fn collect_blocks(&mut self, scope: usize, node: Node<'a>) {
        if node.kind() == "block_construct" {
            self.add_scope(self.scopes[scope].unit, node, Some(scope), false, false);
        } else if !PROGRAM_UNITS.contains(&node.kind()) {
            for item in children(node) {
                self.collect_blocks(scope, item);
            }
        }
    }

    fn modify_variables(&mut self, scope: usize, node: Node<'a>) {
        let unit = self.unit(scope);
        let attrs = attributes(unit, node);
        let procedure = attrs.contains("external") || attrs.contains("intrinsic");
        for item in fields(node, "declarator") {
            let name = unit.name(declaration_name(item));
            if let Some(&id) = self.scopes[scope].symbols.get(&name) {
                self.symbols[id].attrs.extend(attrs.iter().cloned());
                if procedure {
                    self.symbols[id].kind = "procedure";
                }
                if let Some(shape) = declared_shape(unit, item, node) {
                    self.symbols[id].rank = shape_rank(shape);
                }
            } else if procedure {
                let mut symbol = self.symbol(scope, name.clone(), "procedure", item);
                symbol.is_public = self.scopes[scope].public(&name, &attrs);
                symbol.attrs = attrs.clone();
                self.register(symbol, None);
            }
        }
    }

    fn variables(
        &mut self,
        scope: usize,
        node: Node<'a>,
        owner: Option<usize>,
        components_private: bool,
    ) {
        let unit = self.unit(scope);
        let dtype = node.child_by_field_name("type");
        let attrs = attributes(unit, node);
        let kind =
            if dtype.is_some_and(|item| item.kind() == "procedure") || attrs.contains("external") {
                "procedure"
            } else if attrs.contains("parameter") {
                "constant"
            } else {
                "variable"
            };
        for item in fields(node, "declarator") {
            let name_node = declaration_name(item);
            let name = unit.name(name_node);
            let mut symbol = self.symbol(scope, name.clone(), kind, name_node);
            symbol.attrs = attrs.clone();
            symbol.is_public = if owner.is_some() {
                !attrs.contains("private") && (attrs.contains("public") || !components_private)
            } else {
                self.scopes[scope].public(&name, &attrs)
            };
            symbol.type_name = declared_type(unit, dtype);
            symbol.intrinsic = intrinsic_type(dtype).map(str::to_owned);
            if matches!(kind, "variable" | "constant") {
                symbol.rank = declared_shape(unit, item, node).map_or(Some(0), shape_rank);
            }
            symbol.interface_name = dtype
                .and_then(|part| child(part, "procedure_interface"))
                .map(|part| unit.name(part));
            self.register(symbol, owner);
        }
    }

    fn derived_type(&mut self, scope: usize, node: Node<'a>) {
        let Some(header) = child(node, "derived_type_statement") else {
            return;
        };
        let Some(name_node) = child(header, "type_name") else {
            return;
        };
        let unit = self.unit(scope);
        let name = unit.name(name_node);
        let mut symbol = self.symbol(scope, name.clone(), "type", name_node);
        symbol.attrs = attributes(unit, header);
        symbol.is_public = self.scopes[scope].public(&name, &symbol.attrs);
        let base = header.child_by_field_name("base");
        symbol.base_name = base
            .and_then(|item| item.named_child(0))
            .map(|item| unit.name(item));
        let base_name = symbol.base_name.clone();
        let type_id = self.register(symbol, None);
        if !self.types.contains(&type_id) && self.symbols[type_id].kind == "type" {
            self.types.push(type_id);
        }
        let components_private = child(node, "private_statement").is_some();
        for item in children(node) {
            match item.kind() {
                "variable_declaration" => {
                    self.variables(scope, item, Some(type_id), components_private)
                }
                "derived_type_procedures" => {
                    let private = child(item, "private_statement").is_some();
                    for binding in children(item) {
                        self.binding(type_id, binding, private);
                    }
                }
                _ => {}
            }
        }
        if let (Some(base), Some(name)) = (base, base_name)
            && !self.symbols[type_id].members.contains_key(&name)
        {
            let mut parent = self.symbol(scope, name.clone(), "parent_component", base);
            parent.type_name = Some(name);
            parent.is_public = !components_private;
            parent.rank = Some(0);
            self.register(parent, Some(type_id));
        }
    }

    fn binding(&mut self, owner: usize, node: Node<'a>, private: bool) {
        let kind = match node.kind() {
            "procedure_statement" => "binding",
            "generic_statement" => "generic",
            "final_statement" => "final",
            _ => return,
        };
        let scope = self.symbols[owner].scope;
        let unit = self.unit(scope);
        let attrs = attributes(unit, node);
        let interface = child(node, "procedure_interface");
        for item in fields(node, "declarator") {
            let boxed = matches!(item.kind(), "binding" | "binding_list");
            let name_node = if boxed {
                let Some(part) = child(item, "binding_name").and_then(|name| name.named_child(0))
                else {
                    continue;
                };
                part
            } else {
                item
            };
            let name = unit.name(name_node);
            let mut symbol = self.symbol(scope, name.clone(), kind, name_node);
            symbol.is_public = !attrs.contains("private") && (attrs.contains("public") || !private);
            symbol.attrs = attrs.clone();
            let mut specifics = if boxed {
                children(item)
                    .into_iter()
                    .filter(|part| part.kind() == "method_name")
                    .map(|part| unit.name(part))
                    .collect()
            } else {
                Vec::new()
            };
            match kind {
                "binding" => {
                    symbol.target = specifics
                        .pop()
                        .or_else(|| {
                            if attrs.contains("deferred") {
                                interface.map(|part| unit.name(part))
                            } else {
                                None
                            }
                        })
                        .or(Some(name));
                }
                "final" => symbol.target = Some(name),
                _ => symbol.specifics = specifics,
            }
            self.register(symbol, Some(owner));
        }
    }

    fn interface(&mut self, scope: usize, node: Node<'a>) {
        let Some(header) = child(node, "interface_statement") else {
            return;
        };
        let unit = self.unit(scope);
        let abstract_ = child(header, "abstract_specifier").is_some();
        let named = children(header)
            .into_iter()
            .find(|part| NAME_NODES.contains(&part.kind()));
        let mut specifics = Vec::new();
        for item in children(node) {
            if matches!(item.kind(), "function" | "subroutine") {
                let body =
                    self.add_scope(self.scopes[scope].unit, item, Some(scope), true, abstract_);
                specifics.push(self.scopes[body].name.clone());
            } else if item.kind() == "procedure_statement" {
                specifics.extend(
                    fields(item, "declarator")
                        .into_iter()
                        .map(|part| unit.name(declaration_name(part))),
                );
            }
        }
        if let Some(node) = named {
            let name = unit.name(node);
            if let Some(&id) = self.scopes[scope].generics.get(&name) {
                self.symbols[id].specifics.extend(specifics);
            } else {
                let mut symbol = self.symbol(scope, name.clone(), "generic", node);
                symbol.is_public = self.scopes[scope].public(&name, &symbol.attrs);
                symbol.specifics = specifics;
                let id = self.push_symbol(symbol);
                self.scopes[scope].generics.insert(name, id);
            }
        }
    }

    fn link_implementations(&mut self) {
        for scope in 0..self.scopes.len() {
            let item = &self.scopes[scope];
            if item.kind != "submodule" {
                continue;
            }
            let ancestor = item.ancestor.as_deref().unwrap_or("");
            let host = match &item.parent_name {
                Some(parent) => self.submodules.get(&format!("{ancestor}:{parent}")),
                None => self.modules.get(ancestor),
            }
            .copied();
            if let Some(host) = host {
                self.hosts.insert(scope, host);
            } else {
                let name = item.parent_name.as_deref().unwrap_or(ancestor).to_owned();
                self.issue(
                    "SEM004",
                    scope,
                    item.header,
                    &name,
                    format!("Unresolved submodule host '{name}'"),
                );
            }
        }
        for scope in 0..self.scopes.len() {
            let item = &self.scopes[scope];
            let Some(parent) = item.parent else { continue };
            if item.interface_body {
                if !item.attrs.contains("module")
                    && !item.attrs.contains("bind(c)")
                    && item
                        .symbol
                        .is_some_and(|id| self.symbols[id].kind != "abstract_interface")
                    && let Some(&external) = self.external_bodies.get(&item.name)
                {
                    self.implementations.insert(scope, external);
                }
                continue;
            }
            if !matches!(self.scopes[parent].kind, "module" | "submodule")
                || (item.kind != "module_procedure" && !item.attrs.contains("module"))
            {
                continue;
            }
            let name = item.name.clone();
            let header = item.header;
            let in_submodule = self.scopes[parent].kind == "submodule";
            let mut host = if in_submodule {
                self.host(parent)
            } else {
                Some(parent)
            };
            let mut seen = HashSet::new();
            let mut declaration = None;
            while let Some(candidate) = host.filter(|id| seen.insert(*id)) {
                let declared = self.scopes[candidate]
                    .symbols
                    .get(&name)
                    .and_then(|id| self.symbols[*id].procedure);
                if declared.is_some_and(|id| {
                    self.scopes[id].interface_body && self.scopes[id].attrs.contains("module")
                }) {
                    declaration = declared;
                    break;
                }
                host = self.host(candidate);
            }
            let Some(declaration) = declaration else {
                if in_submodule || item.kind == "module_procedure" {
                    self.issue(
                        "SEM004",
                        scope,
                        header,
                        &name,
                        "Separate procedure has no ancestor interface".into(),
                    );
                }
                continue;
            };
            if self.implementations.contains_key(&declaration) {
                self.issue(
                    "SEM004",
                    scope,
                    header,
                    &name,
                    "Separate procedure has multiple implementations".into(),
                );
                continue;
            }
            self.implementations.insert(declaration, scope);
            self.scopes[scope].declaration = Some(declaration);
            let parameters = self.scopes[declaration].parameters.clone();
            if self.scopes[scope].kind == "module_procedure" {
                let result = self.scopes[declaration].result.clone();
                self.scopes[scope].parameters = parameters.clone();
                self.scopes[scope].result = result.clone();
                for argument in parameters.into_iter().chain(result) {
                    if let Some(&id) = self.scopes[declaration].symbols.get(&argument) {
                        self.scopes[scope].symbols.entry(argument).or_insert(id);
                    }
                }
            } else if self.scopes[scope].parameters != parameters {
                self.issue(
                    "SEM004",
                    scope,
                    header,
                    &name,
                    "Separate procedure arguments do not match its ancestor interface".into(),
                );
            }
        }
        for scope in 0..self.scopes.len() {
            let item = &self.scopes[scope];
            if item.interface_body
                && item.attrs.contains("module")
                && !self.implementations.contains_key(&scope)
            {
                let name = item.name.clone();
                self.issue(
                    "SEM004",
                    scope,
                    item.header,
                    &name,
                    format!("Separate procedure '{name}' has no implementation"),
                );
            }
        }
    }

    fn local(&self, scope: usize, name: &str, kinds: &[&str]) -> Option<usize> {
        let scope = &self.scopes[scope];
        scope
            .generics
            .get(name)
            .into_iter()
            .chain(scope.symbols.get(name))
            .copied()
            .find(|id| kinds.is_empty() || kinds.contains(&self.symbols[*id].kind))
    }

    fn candidate_from_use(
        &self,
        use_: &Use,
        name: &str,
        kinds: &[&str],
        active: &mut HashSet<(usize, String, bool)>,
    ) -> Vec<Candidate> {
        let Some(remote) = use_.remote(name) else {
            return Vec::new();
        };
        if !use_.intrinsic
            && let Some(&module) = self.modules.get(&use_.module)
        {
            return self.candidates(module, remote, kinds, true, active);
        }
        let Some(external) = self.externals.get(&use_.module) else {
            return Vec::new();
        };
        if (use_.intrinsic && !external.intrinsic)
            || (use_.non_intrinsic && external.intrinsic)
            || (!external.compiler && !external.symbols.contains(remote))
        {
            return Vec::new();
        }
        vec![Candidate::External(use_.module.clone(), remote.into())]
    }

    fn opaque(&self, candidate: &Candidate) -> bool {
        matches!(candidate, Candidate::External(module, _) if self.externals.get(module).is_some_and(|entry| entry.compiler))
    }

    fn candidates(
        &self,
        scope: usize,
        name: &str,
        kinds: &[&str],
        exported: bool,
        active: &mut HashSet<(usize, String, bool)>,
    ) -> Vec<Candidate> {
        let key = (scope, name.to_owned(), exported);
        if !active.insert(key.clone()) {
            return Vec::new();
        }
        let result = self.search(scope, name, kinds, exported, active);
        active.remove(&key);
        result
    }

    fn search(
        &self,
        scope: usize,
        name: &str,
        kinds: &[&str],
        exported: bool,
        active: &mut HashSet<(usize, String, bool)>,
    ) -> Vec<Candidate> {
        let entry = &self.scopes[scope];
        if let Some(id) = self.local(scope, name, kinds) {
            return if !exported || self.symbols[id].is_public {
                vec![Candidate::Local(id)]
            } else {
                Vec::new()
            };
        }
        if exported && !entry.public(name, &BTreeSet::new()) {
            return Vec::new();
        }
        let mut certain = BTreeSet::new();
        let mut opaque = BTreeSet::new();
        for use_ in &entry.uses {
            for item in self.candidate_from_use(use_, name, kinds, active) {
                if self.opaque(&item)
                    && !use_.imports.contains_key(name)
                    && !(exported && entry.explicit_access.contains_key(name))
                {
                    opaque.insert(item);
                } else {
                    certain.insert(item);
                }
            }
        }
        if !certain.is_empty() {
            return certain.into_iter().collect();
        }
        if !exported
            && entry
                .imports
                .as_ref()
                .is_none_or(|imports| imports.contains(name))
            && let Some(host) = self.host(scope)
        {
            let found = self.candidates(host, name, kinds, false, active);
            if !found.is_empty() {
                return found;
            }
        }
        opaque.into_iter().collect()
    }

    fn materialize(&mut self, scope: usize, node: Node<'a>, candidate: Candidate) -> usize {
        match candidate {
            Candidate::Local(id) => id,
            Candidate::External(module, name) => {
                let key = (module.clone(), name.clone());
                if let Some(&id) = self.external_symbols.get(&key) {
                    return id;
                }
                let kind = if module == "<intrinsic>" {
                    "intrinsic"
                } else {
                    "external"
                };
                let mut symbol = self.symbol(scope, name.clone(), kind, node);
                symbol.is_public = true;
                symbol.qualified = format!("{module}::{name}");
                symbol.external = Some(module);
                let id = self.push_symbol(symbol);
                self.external_symbols.insert(key, id);
                id
            }
        }
    }

    fn named_use(&self, scope: usize, name: &str, visited: &mut HashSet<(usize, String)>) -> bool {
        if !visited.insert((scope, name.into())) {
            return false;
        }
        let entry = &self.scopes[scope];
        for use_ in &entry.uses {
            if use_.imports.contains_key(name) {
                return true;
            }
            if let (Some(&module), Some(remote)) =
                (self.modules.get(&use_.module), use_.remote(name))
                && self.scopes[module].public(remote, &BTreeSet::new())
                && (self.scopes[module].explicit_access.contains_key(remote)
                    || self.named_use(module, remote, visited))
            {
                return true;
            }
        }
        entry
            .imports
            .as_ref()
            .is_none_or(|imports| imports.contains(name))
            && self
                .host(scope)
                .is_some_and(|host| self.named_use(host, name, visited))
    }

    pub fn lookup(&self, scope: usize, name: &str) -> Option<usize> {
        self.lookup_kind(scope, name, &[])
    }

    pub fn lookup_kind(&self, scope: usize, name: &str, kinds: &[&str]) -> Option<usize> {
        let found = self.candidates(
            scope,
            &name.to_ascii_lowercase(),
            kinds,
            false,
            &mut HashSet::new(),
        );
        match found.as_slice() {
            [Candidate::Local(id)] => Some(*id),
            [Candidate::External(module, name)] => self
                .external_symbols
                .get(&(module.clone(), name.clone()))
                .copied(),
            _ => None,
        }
    }

    pub fn resolve_parent(&self, owner: usize) -> Option<usize> {
        let symbol = &self.symbols[owner];
        symbol
            .base_name
            .as_deref()
            .and_then(|name| self.lookup_kind(symbol.scope, name, TYPES))
    }

    pub fn member(&self, mut owner: usize, name: &str) -> Option<usize> {
        let name = name.to_ascii_lowercase();
        let mut visited = HashSet::new();
        while visited.insert(owner) {
            if let Some(&member) = self.symbols[owner].members.get(&name) {
                return Some(member);
            }
            owner = self.resolve_parent(owner)?;
        }
        None
    }

    pub fn binding_target(&self, binding: usize) -> Option<usize> {
        let symbol = &self.symbols[binding];
        let kinds = if symbol.attrs.contains("deferred") {
            INTERFACES
        } else {
            PROCEDURES
        };
        symbol
            .target
            .as_ref()
            .and_then(|target| self.lookup_kind(symbol.scope, target, kinds))
    }

    pub fn implementation(&self, procedure: usize) -> Option<usize> {
        let scope = self.symbols[procedure].procedure?;
        if self.scopes[scope].interface_body {
            self.implementations.get(&scope).copied()
        } else {
            Some(scope)
        }
    }

    pub fn expression(&self, unit: usize, node: Node<'_>) -> Option<usize> {
        self.expression_symbols.get(&(unit, node.id())).copied()
    }

    fn result_of(&self, procedure: usize) -> Option<usize> {
        let scope = &self.scopes[self.symbols[procedure].procedure?];
        scope.symbols.get(scope.result.as_ref()?).copied()
    }

    pub fn type_of(&self, symbol: usize) -> Option<usize> {
        self.type_of_inner(symbol, &mut HashSet::new())
    }

    fn type_of_inner(&self, symbol: usize, visited: &mut HashSet<usize>) -> Option<usize> {
        if !visited.insert(symbol) {
            return None;
        }
        let entry = &self.symbols[symbol];
        let result = match entry.kind {
            "type" => Some(symbol),
            "binding" => self
                .binding_target(symbol)
                .and_then(|target| self.type_of_inner(target, visited)),
            "generic" => {
                let mut types = self
                    .generic_specifics(symbol, None)
                    .into_iter()
                    .map(|id| self.type_of_inner(id, visited));
                let first = types.next().flatten();
                first.filter(|id| types.all(|item| item == Some(*id)))
            }
            _ => self
                .result_of(symbol)
                .and_then(|id| self.type_of_inner(id, visited))
                .or_else(|| {
                    entry
                        .type_name
                        .as_deref()
                        .filter(|name| *name != "*")
                        .and_then(|name| self.lookup_kind(entry.scope, name, TYPES))
                }),
        };
        visited.remove(&symbol);
        result
    }

    pub fn expression_category(&self, unit: usize, node: Node<'_>) -> Option<&str> {
        match node.kind() {
            "boolean_literal" | "logical_expression" | "relational_expression" => {
                return Some("logical");
            }
            "string_literal" => return Some("character"),
            "complex_literal" => return Some("complex"),
            "number_literal" => {
                let spell = self.units[unit].text(node).split('_').next().unwrap_or("");
                if !spell.is_empty() && spell.bytes().all(|byte| byte.is_ascii_digit()) {
                    return Some("integer");
                }
                if spell
                    .as_bytes()
                    .first()
                    .is_some_and(|byte| byte.is_ascii_digit() || *byte == b'.')
                {
                    return Some("real");
                }
            }
            "unary_expression" | "parenthesized_expression" | "array_literal" => {
                if let Some(inner) = children(node)
                    .into_iter()
                    .find(|part| part.kind() != "comment")
                {
                    return self.expression_category(unit, inner);
                }
            }
            _ => {}
        }
        let mut symbol = self.expression(unit, node)?;
        if self.symbols[symbol].kind == "binding" {
            symbol = self.binding_target(symbol)?
        }
        self.symbols[self.result_of(symbol).unwrap_or(symbol)]
            .intrinsic
            .as_deref()
    }

    fn expression_rank(&self, unit: usize, node: Node<'_>) -> Option<usize> {
        match node.kind() {
            "number_literal" | "boolean_literal" | "string_literal" | "complex_literal" => {
                return Some(0);
            }
            "array_literal" => return Some(1),
            "unary_expression" | "parenthesized_expression" => {
                return self.expression_rank(unit, node.named_child(0)?);
            }
            "math_expression" | "relational_expression" | "logical_expression" => {
                let left = self.expression_rank(unit, node.child_by_field_name("left")?)?;
                let right = self.expression_rank(unit, node.child_by_field_name("right")?)?;
                return Some(left.max(right));
            }
            _ => {}
        }
        let mut id = self.expression(unit, node)?;
        let symbol = &self.symbols[id];
        if symbol.kind == "type" && node.kind() == "call_expression" {
            return Some(0);
        }
        if symbol.kind == "binding" {
            id = self.binding_target(id)?;
        }
        let symbol = &self.symbols[id];
        if symbol.attrs.contains("elemental") {
            return None;
        }
        let rank = self.symbols[self.result_of(id).unwrap_or(id)].rank?;
        if node.kind() == "call_expression"
            && rank != 0
            && matches!(symbol.kind, "variable" | "constant")
        {
            let subscripts = child(node, "argument_list")?;
            return children(subscripts)
                .into_iter()
                .filter(|part| part.kind() != "comment")
                .try_fold(0, |rank, subscript| {
                    if subscript.kind() == "extent_specifier" {
                        Some(rank + 1)
                    } else {
                        self.expression_rank(unit, subscript)
                            .map(|n| rank + usize::from(n != 0))
                    }
                });
        }
        if node.kind() == "derived_type_member_expression" {
            let receiver_rank = self.expression_rank(unit, node.named_child(0)?)?;
            return Some(rank.max(receiver_rank));
        }
        Some(rank)
    }

    fn check_declarations(&mut self) {
        for scope in 0..self.scopes.len() {
            let header = self.scopes[scope].header;
            for use_ in self.scopes[scope].uses.clone() {
                for (local, remote) in &use_.imports {
                    let candidates =
                        self.candidate_from_use(&use_, local, &[], &mut HashSet::new());
                    let location = use_.imported_at[local].clone();
                    if candidates.is_empty() {
                        let subject = format!("{}::{local}", self.scopes[scope].qualified);
                        self.diagnostics.push(
                            Diagnostic::new(
                                "SEM002",
                                format!(
                                    "'{remote}' is not an exported symbol of '{}'",
                                    use_.module
                                ),
                                location.clone(),
                            )
                            .with_subject(subject),
                        );
                    }
                    for candidate in candidates {
                        if matches!(candidate, Candidate::External(..)) {
                            let external = self.materialize(scope, header, candidate);
                            self.external_imports
                                .push((scope, external, location.clone()));
                        }
                    }
                }
            }
            for (name, node) in self.scopes[scope].explicit_access_nodes.clone() {
                self.resolve(scope, node, &name, &[], false);
            }
            let generics: Vec<_> = self.scopes[scope].generics.values().copied().collect();
            for id in generics {
                let (name, node) = (self.symbols[id].name.clone(), self.symbols[id].node);
                for specific in self.symbols[id].specifics.clone() {
                    if !self.resolve_declaration(scope, node, &specific, PROCEDURES) {
                        self.issue(
                            "SEM005",
                            scope,
                            node,
                            &specific,
                            format!(
                                "Generic '{name}' has unresolved specific procedure '{specific}'"
                            ),
                        );
                    }
                }
            }
            if let Some(imports) = self.scopes[scope].imports.clone() {
                for name in imports {
                    if !self
                        .host(scope)
                        .is_some_and(|host| self.resolve_declaration(host, header, &name, &[]))
                    {
                        self.issue(
                            "SEM002",
                            scope,
                            header,
                            &name,
                            format!("Unresolved host import '{name}'"),
                        );
                    }
                }
            }
            let arguments = self.scopes[scope]
                .parameters
                .clone()
                .into_iter()
                .chain(self.scopes[scope].result.clone());
            for name in arguments {
                if !self.scopes[scope].symbols.contains_key(&name) {
                    self.issue(
                        "SEM002",
                        scope,
                        header,
                        &name,
                        format!("Argument or function result '{name}' is undeclared"),
                    );
                }
            }
        }
        for owner in self.types.clone() {
            let scope = self.symbols[owner].scope;
            if let Some(name) = self.symbols[owner].base_name.clone() {
                let node = self.symbols[owner].node;
                self.resolve_declaration(scope, node, &name, TYPES);
            }
            let members: Vec<_> = self.symbols[owner].members.values().copied().collect();
            for id in members {
                let kind = self.symbols[id].kind;
                let node = self.symbols[id].node;
                let name = self.symbols[id].name.clone();
                if matches!(kind, "binding" | "final")
                    && let Some(target) = self.symbols[id].target.clone()
                {
                    let kinds = if self.symbols[id].attrs.contains("deferred") {
                        INTERFACES
                    } else {
                        PROCEDURES
                    };
                    if !self.resolve_declaration(scope, node, &target, kinds) {
                        let owner_name = &self.symbols[owner].name;
                        self.issue("SEM005", scope, node, &name,
                            format!("Unresolved implementation/interface '{target}' of '{owner_name}%{name}'"));
                    }
                } else if kind == "generic" {
                    for specific in self.symbols[id].specifics.clone() {
                        if self
                            .member(owner, &specific)
                            .is_none_or(|item| self.symbols[item].kind != "binding")
                        {
                            let owner_name = &self.symbols[owner].name;
                            self.issue("SEM005", scope, node, &specific,
                                format!("Unresolved specific binding '{specific}' of '{owner_name}%{name}'"));
                        }
                    }
                }
            }
        }
    }

    fn resolve_declaration(
        &mut self,
        scope: usize,
        node: Node<'a>,
        name: &str,
        kinds: &[&str],
    ) -> bool {
        let found = self.candidates(scope, name, kinds, false, &mut HashSet::new());
        if let [candidate] = found.as_slice() {
            self.materialize(scope, node, candidate.clone());
            true
        } else {
            false
        }
    }

    fn resolve(
        &mut self,
        scope: usize,
        node: Node<'a>,
        name: &str,
        kinds: &[&str],
        intrinsic: bool,
    ) -> Option<usize> {
        let found = self.candidates(scope, name, kinds, false, &mut HashSet::new());
        let opaque = found.iter().all(|candidate| self.opaque(candidate));
        if intrinsic
            && self.intrinsics.contains(name)
            && opaque
            && (found.is_empty() || !self.named_use(scope, name, &mut HashSet::new()))
        {
            return Some(self.materialize(
                scope,
                node,
                Candidate::External("<intrinsic>".into(), name.into()),
            ));
        }
        if name.starts_with("mpi_") && self.externals.contains_key("mpi")
            && !found.is_empty() && opaque
            && found.iter().all(|candidate| matches!(candidate, Candidate::External(module, _) if module.starts_with("petsc")))
        {
            return Some(self.materialize(scope, node, Candidate::External("mpi".into(), name.into())));
        }
        match found.as_slice() {
            [candidate] => Some(self.materialize(scope, node, candidate.clone())),
            [] => {
                self.issue(
                    "SEM002",
                    scope,
                    node,
                    name,
                    format!("Unresolved symbol '{name}'"),
                );
                None
            }
            _ if opaque => None,
            _ => {
                let origins = found
                    .iter()
                    .map(|candidate| match candidate {
                        Candidate::Local(id) => self.symbols[*id].qualified.clone(),
                        Candidate::External(module, name) => format!("{module}::{name}"),
                    })
                    .collect::<Vec<_>>()
                    .join(", ");
                self.issue(
                    "SEM003",
                    scope,
                    node,
                    name,
                    format!("Ambiguous name '{name}' ({origins})"),
                );
                None
            }
        }
    }

    fn generic_specifics(&self, generic: usize, dispatch: Option<usize>) -> Vec<usize> {
        let generic = &self.symbols[generic];
        generic
            .specifics
            .iter()
            .filter_map(|name| {
                if let Some(owner) = dispatch.or(generic.owner) {
                    self.member(owner, name)
                } else {
                    self.lookup_kind(generic.scope, name, PROCEDURES)
                }
            })
            .filter(|id| self.symbols[*id].kind != "generic")
            .collect()
    }

    fn select_call(&self, symbol: usize, actuals: &[Argument], dispatch: Option<usize>) -> usize {
        if self.symbols[symbol].kind != "generic" {
            return symbol;
        }
        let specifics = self.generic_specifics(symbol, dispatch);
        if specifics.len() != self.symbols[symbol].specifics.len() {
            return symbol;
        }
        let candidates: Vec<_> = specifics
            .into_iter()
            .filter(|id| self.matches_arguments(*id, actuals))
            .collect();
        if let [specific] = candidates.as_slice() {
            *specific
        } else {
            symbol
        }
    }

    fn matches_arguments(&self, specific: usize, actuals: &[Argument]) -> bool {
        let binding = &self.symbols[specific];
        let target = if binding.kind == "binding" {
            self.binding_target(specific)
        } else {
            Some(specific)
        };
        let Some(scope) = target.and_then(|id| self.symbols[id].procedure) else {
            return true;
        };
        let signature = &self.scopes[scope];
        let mut parameters: Vec<_> = signature.parameters.iter().map(String::as_str).collect();
        if binding.kind == "binding" && !parameters.is_empty() && !binding.attrs.contains("nopass")
        {
            let passed = binding
                .attrs
                .iter()
                .find_map(|attr| attr.strip_prefix("pass(")?.strip_suffix(')'))
                .unwrap_or(parameters[0]);
            parameters.retain(|name| *name != passed);
        }
        let mut matched = BTreeSet::new();
        for (index, actual) in actuals.iter().enumerate() {
            let Some(name) = actual
                .keyword
                .as_deref()
                .or_else(|| parameters.get(index).copied())
            else {
                return false;
            };
            if !parameters.contains(&name) || !matched.insert(name) {
                return false;
            }
            if let Some(&id) = signature.symbols.get(name) {
                let dummy = &self.symbols[id];
                if let (Some(expected), Some(observed)) = (dummy.rank, actual.rank)
                    && expected != observed
                    && !(expected == 0 && signature.attrs.contains("elemental"))
                {
                    return false;
                }
                if let (Some(intrinsic), Some(category)) = (&dummy.intrinsic, &actual.category)
                    && intrinsic != category
                {
                    return false;
                }
                if let (Some(actual), Some(required)) = (actual.derived, self.type_of(id))
                    && actual != required
                {
                    let mut current = Some(actual);
                    let mut seen = HashSet::new();
                    let mut found = false;
                    while let Some(owner) = current.filter(|id| seen.insert(*id)) {
                        if owner == required {
                            found = true;
                            break;
                        }
                        current = self.resolve_parent(owner);
                    }
                    if !found {
                        return false;
                    }
                }
            }
        }
        parameters.into_iter().all(|name| {
            matched.contains(name)
                || signature
                    .symbols
                    .get(name)
                    .is_none_or(|id| self.symbols[*id].attrs.contains("optional"))
        })
    }

    fn record(
        &mut self,
        scope: usize,
        node: Node<'a>,
        symbol: usize,
        receiver: Option<usize>,
        root: Option<usize>,
        call: bool,
    ) {
        self.references.push(Reference {
            scope,
            node,
            symbol,
            receiver,
            root_receiver: root,
            call,
        });
        self.expression_symbols
            .insert((self.scopes[scope].unit, node.id()), symbol);
    }

    fn scan_body(&mut self, scope: usize, aliases: &Aliases) {
        let header = self.scopes[scope].header;
        for item in children(header)
            .into_iter()
            .filter(|item| matches!(item.kind(), "intrinsic_type" | "derived_type"))
        {
            self.scan(scope, item, aliases, None, false);
        }
        for item in children(self.scopes[scope].node) {
            if item != header && !PROGRAM_UNITS.contains(&item.kind()) {
                self.scan(scope, item, aliases, None, false);
            }
        }
    }

    fn scan(
        &mut self,
        scope: usize,
        node: Node<'a>,
        aliases: &Aliases,
        actuals: Option<&[Argument]>,
        subroutine: bool,
    ) -> Option<usize> {
        let unit_index = self.scopes[scope].unit;
        let unit = &self.units[unit_index];
        match node.kind() {
            "identifier" => {
                let name = unit.name(node);
                let symbol = if let Some(&alias) = aliases.get(&name) {
                    alias
                } else {
                    self.resolve(
                        scope,
                        node,
                        &name,
                        if subroutine { CALLABLE } else { &[] },
                        actuals.is_some(),
                    )
                }?;
                let symbol = actuals.map_or(symbol, |args| self.select_call(symbol, args, None));
                let kind = self.symbols[symbol].kind;
                let call = actuals.is_some()
                    && (CALLABLE.contains(&kind) || matches!(kind, "external" | "intrinsic"));
                self.record(scope, node, symbol, None, None, call);
                Some(symbol)
            }
            "derived_type_member_expression" => self.scan_member(scope, node, aliases, actuals),
            "call_expression" | "subroutine_call" => {
                let named = children(node);
                let target = *named.iter().find(|item| item.kind() != "comment")?;
                let mut args = Vec::new();
                if let Some(list) = child(node, "argument_list") {
                    for item in children(list)
                        .into_iter()
                        .filter(|part| part.kind() != "comment")
                    {
                        let key = (item.kind() == "keyword_argument")
                            .then(|| item.child_by_field_name("name"))
                            .flatten();
                        let value = if key.is_some() {
                            item.child_by_field_name("value")
                        } else {
                            Some(item)
                        };
                        if let Some(value) = value {
                            let observed = self.scan(scope, value, aliases, None, false);
                            let category = self
                                .expression_category(unit_index, value)
                                .map(str::to_owned);
                            let derived = observed.and_then(|id| self.type_of(id));
                            args.push(Argument {
                                keyword: key.map(|key| unit.name(key)),
                                category,
                                derived,
                                rank: self.expression_rank(unit_index, value),
                            });
                        }
                    }
                }
                let symbol = self.scan(
                    scope,
                    target,
                    aliases,
                    Some(&args),
                    node.kind() == "subroutine_call",
                )?;
                self.expression_symbols
                    .insert((unit_index, node.id()), symbol);
                if let Some(&root) = self.receiver_roots.get(&(unit_index, target.id())) {
                    self.receiver_roots.insert((unit_index, node.id()), root);
                }
                Some(symbol)
            }
            "keyword_argument" => node
                .child_by_field_name("value")
                .and_then(|value| self.scan(scope, value, aliases, None, false)),
            "variable_declaration" | "variable_modification" => {
                let declared = fields(node, "declarator");
                for item in children(node) {
                    if declared.contains(&item) {
                        self.scan_declarator(scope, item, aliases)
                    } else {
                        self.scan(scope, item, aliases, None, false);
                    }
                }
                None
            }
            "derived_type" | "procedure_interface" | "base_type_specifier" => {
                self.scan_type(scope, node, aliases)
            }
            "derived_type_definition" => {
                if let Some(header) = child(node, "derived_type_statement") {
                    if let Some(base) = header.child_by_field_name("base") {
                        self.scan_type(scope, base, aliases);
                    }
                    let mut type_aliases = aliases.clone();
                    if let Some(name) = child(header, "type_name").map(|item| unit.name(item))
                        && let Some(&owner) = self.scopes[scope].symbols.get(&name)
                    {
                        for &id in self.symbols[owner].members.values() {
                            let item = &self.symbols[id];
                            if item.attrs.contains("kind") || item.attrs.contains("len") {
                                type_aliases.insert(item.name.clone(), Some(id));
                            }
                        }
                    }
                    for item in children(node)
                        .into_iter()
                        .filter(|part| part.kind() == "variable_declaration")
                    {
                        self.scan(scope, item, &type_aliases, None, false);
                    }
                }
                None
            }
            "associate_statement"
            | "select_type_statement"
            | "select_rank_statement"
            | "implied_do_loop_expression" => {
                self.scan_associations(scope, node, aliases);
                None
            }
            "block_construct" => {
                if let Some(&block) = self.scope_by_node.get(&(unit_index, node.id())) {
                    let visible = aliases
                        .iter()
                        .filter(|(name, _)| !self.scopes[block].symbols.contains_key(*name))
                        .map(|(name, value)| (name.clone(), *value))
                        .collect();
                    self.scan_body(block, &visible);
                }
                None
            }
            kind if SKIP.contains(&kind)
                || kind.starts_with("end_")
                || PROGRAM_UNITS.contains(&kind) =>
            {
                None
            }
            _ => {
                for item in children(node) {
                    self.scan(scope, item, aliases, None, false);
                }
                None
            }
        }
    }

    fn scan_member(
        &mut self,
        scope: usize,
        node: Node<'a>,
        aliases: &Aliases,
        actuals: Option<&[Argument]>,
    ) -> Option<usize> {
        let unit_index = self.scopes[scope].unit;
        let unit = &self.units[unit_index];
        let names: Vec<_> = children(node)
            .into_iter()
            .filter(|part| part.kind() != "comment")
            .collect();
        let [receiver, member] = names.as_slice() else {
            return None;
        };
        let receiver_id = self.scan(scope, *receiver, aliases, None, false)?;
        let root = self
            .receiver_roots
            .get(&(unit_index, receiver.id()))
            .copied()
            .unwrap_or(receiver_id);
        self.receiver_roots.insert((unit_index, node.id()), root);
        let owner = self.type_of(receiver_id)?;
        if self.symbols[owner].external.is_some() {
            return None;
        }
        let name = unit.name(*member);
        let Some(id) = self.member(owner, &name) else {
            let mut current = Some(owner);
            let mut visited = HashSet::new();
            while let Some(parent) = current.filter(|id| visited.insert(*id)) {
                if self.symbols[parent].external.is_some() {
                    return None;
                }
                current = self.resolve_parent(parent);
            }
            let owner = &self.symbols[owner].qualified;
            self.issue(
                "SEM002",
                scope,
                *member,
                &name,
                format!("No member '{name}' of '{owner}'"),
            );
            return None;
        };
        let symbol = &self.symbols[id];
        if self.owner_scope(scope) != self.owner_scope(symbol.scope) {
            let visibility = if symbol.kind == "parent_component" {
                self.type_of(id)
            } else {
                Some(id)
            };
            if visibility.is_some_and(|id| !self.symbols[id].is_public) {
                let member_name = &symbol.qualified;
                self.issue(
                    "SEM003",
                    scope,
                    *member,
                    &name,
                    format!("Private member '{member_name}' is inaccessible"),
                );
            }
        }
        let id = actuals.map_or(id, |args| self.select_call(id, args, Some(owner)));
        self.expression_symbols.insert((unit_index, node.id()), id);
        self.record(
            scope,
            *member,
            id,
            Some(receiver_id),
            Some(root),
            actuals.is_some() && CALLABLE.contains(&self.symbols[id].kind),
        );
        Some(id)
    }

    fn scan_declarator(&mut self, scope: usize, node: Node<'a>, aliases: &Aliases) {
        let name = declaration_name(node);
        let unit = self.scopes[scope].unit;
        if let Some(id) = self.expression(unit, name) {
            self.expression_symbols.insert((unit, node.id()), id);
        }
        if name == node {
            return;
        }
        for item in children(node) {
            if item == name {
                continue;
            }
            if item.start_byte() <= name.start_byte() && name.start_byte() < item.end_byte() {
                self.scan_declarator(scope, item, aliases);
            } else {
                self.scan(scope, item, aliases, None, false);
            }
        }
    }

    fn scan_type(&mut self, scope: usize, node: Node<'a>, aliases: &Aliases) -> Option<usize> {
        let (name_node, kinds) = match node.kind() {
            "procedure_interface" => (Some(node), INTERFACES),
            "derived_type" => (
                node.child_by_field_name("name")
                    .filter(|name| name.kind() == "type_name"),
                TYPES,
            ),
            _ => (node.named_child(0), TYPES),
        };
        let mut found = None;
        if let Some(name_node) = name_node {
            let name = self.unit(scope).name(name_node);
            if node.kind() == "procedure_interface"
                && let Some(category) = procedure_type_keyword(node)
                && let Some(declaration) = node
                    .parent()
                    .filter(|item| item.kind() == "procedure")
                    .and_then(|item| item.parent())
                && self
                    .candidates(scope, &name, kinds, false, &mut HashSet::new())
                    .is_empty()
            {
                for item in fields(declaration, "declarator") {
                    if let Some(id) =
                        self.expression(self.scopes[scope].unit, declaration_name(item))
                    {
                        self.symbols[id].interface_name = None;
                        self.symbols[id].intrinsic = Some(category.into());
                    }
                }
                return None;
            }
            if name != "*" && !name.is_empty() {
                found = self.resolve(scope, name_node, &name, kinds, false);
                if let Some(symbol) = found {
                    self.record(scope, name_node, symbol, None, None, false);
                }
            }
        }
        for item in children(node)
            .into_iter()
            .filter(|part| Some(*part) != name_node)
        {
            self.scan(scope, item, aliases, None, false);
        }
        found
    }

    fn scan_associations(&mut self, scope: usize, node: Node<'a>, aliases: &Aliases) {
        let unit = self.unit(scope);
        let mut local = aliases.clone();
        match node.kind() {
            "implied_do_loop_expression" => {
                let control = child(node, "loop_control_expression");
                if let Some(control) = control {
                    if let Some(index) = control.named_child(0) {
                        local.insert(unit.name(index), None);
                    }
                    for item in children(control).into_iter().skip(1) {
                        self.scan(scope, item, &local, None, false);
                    }
                }
                for item in children(node)
                    .into_iter()
                    .filter(|part| Some(*part) != control)
                {
                    self.scan(scope, item, &local, None, false);
                }
            }
            "associate_statement" => {
                let list = child(node, "association_list");
                if let Some(list) = list {
                    for item in children(list) {
                        if let (Some(name), Some(value)) = (
                            item.child_by_field_name("name"),
                            item.child_by_field_name("selector"),
                        ) {
                            local.insert(
                                unit.name(name),
                                self.scan(scope, value, aliases, None, false),
                            );
                        }
                    }
                }
                for item in children(node)
                    .into_iter()
                    .filter(|part| Some(*part) != list)
                {
                    self.scan(scope, item, &local, None, false);
                }
            }
            _ => {
                let selector = child(node, "selector");
                let mut alias = None;
                if let Some(mut value) = selector.and_then(|item| item.named_child(0)) {
                    if value.kind() == "pointer_association_statement" {
                        let parts: Vec<_> = children(value)
                            .into_iter()
                            .filter(|item| item.kind() != "comment")
                            .collect();
                        if let [name, target] = parts.as_slice() {
                            alias = Some(unit.name(*name));
                            value = *target;
                        }
                    } else if value.kind() == "identifier" {
                        alias = Some(unit.name(value));
                    }
                    let target = self.scan(scope, value, aliases, None, false);
                    if let Some(name) = &alias {
                        local.insert(name.clone(), target);
                    }
                }
                for branch in children(node)
                    .into_iter()
                    .filter(|part| Some(*part) != selector)
                {
                    let type_guard = (branch.kind() == "type_statement")
                        .then(|| branch.child_by_field_name("type"))
                        .flatten()
                        .filter(|part| part.kind() == "identifier");
                    if let (Some(guard), Some(name)) = (type_guard, &alias) {
                        let mut scoped = local.clone();
                        let id = self.resolve(scope, guard, &unit.name(guard), TYPES, false);
                        if let Some(id) = id {
                            self.record(scope, guard, id, None, None, false);
                        }
                        scoped.insert(name.clone(), id);
                        for item in children(branch).into_iter().filter(|part| *part != guard) {
                            self.scan(scope, item, &scoped, None, false);
                        }
                    } else {
                        self.scan(scope, branch, &local, None, false);
                    }
                }
            }
        }
    }
}

fn declared_shape<'a>(
    unit: &Unit,
    mut declarator: Node<'a>,
    declaration: Node<'a>,
) -> Option<Node<'a>> {
    while matches!(
        declarator.kind(),
        "sized_declarator" | "coarray_declarator" | "init_declarator" | "pointer_init_declarator"
    ) {
        if let Some(shape) = child(declarator, "size") {
            return Some(shape);
        }
        declarator = declarator
            .child_by_field_name("left")
            .or_else(|| declarator.named_child(0))?;
    }
    children(declaration)
        .into_iter()
        .find(|part| {
            part.kind() == "type_qualifier"
                && part
                    .child(0)
                    .is_some_and(|token| unit.name(token) == "dimension")
        })
        .and_then(|part| child(part, "argument_list"))
}

fn shape_rank(shape: Node<'_>) -> Option<usize> {
    let extents = children(shape);
    if extents.iter().any(|part| part.kind() == "assumed_rank") {
        None
    } else {
        Some(
            extents
                .iter()
                .filter(|part| part.kind() != "comment")
                .count(),
        )
    }
}

pub fn spelling(unit: &Unit, node: Node<'_>) -> String {
    if node.child_count() == 0 {
        return unit.name(node);
    }
    let mut cursor = node.walk();
    node.children(&mut cursor)
        .filter(|part| !matches!(part.kind(), "comment" | "&" | "line_continuation"))
        .map(|part| spelling(unit, part))
        .collect()
}

pub fn attributes(unit: &Unit, node: Node<'_>) -> BTreeSet<String> {
    let mut attrs = BTreeSet::new();
    for item in children(node) {
        match item.kind() {
            "access_specifier"
            | "procedure_attribute"
            | "procedure_qualifier"
            | "type_qualifier" => {
                attrs.insert(spelling(unit, item));
            }
            "abstract_specifier" => {
                attrs.insert("abstract".into());
            }
            "language_binding" if spelling(unit, item).starts_with("bind(c") => {
                attrs.insert("bind(c)".into());
            }
            "procedure" => attrs.extend(attributes(unit, item)),
            _ => {}
        }
    }
    attrs
}

pub fn declaration_name(mut node: Node<'_>) -> Node<'_> {
    while !NAME_NODES.contains(&node.kind()) {
        let Some(next) = node
            .child_by_field_name("left")
            .or_else(|| node.named_child(0))
        else {
            break;
        };
        node = next;
    }
    node
}

pub fn intrinsic_type(node: Option<Node<'_>>) -> Option<&'static str> {
    let node = node?;
    let node = match node.kind() {
        "derived_type" => node.child_by_field_name("name")?,
        "procedure" => child(node, "intrinsic_type")?,
        _ => node,
    };
    if node.kind() != "intrinsic_type" {
        return None;
    }
    let kind = node.child(0)?.kind();
    Some(if kind == "double" { "real" } else { kind })
}

fn declared_type(unit: &Unit, node: Option<Node<'_>>) -> Option<String> {
    let mut node = node?;
    if node.kind() == "procedure" {
        node = child(node, "derived_type")?;
    }
    if node.kind() != "derived_type" {
        return None;
    }
    node.child_by_field_name("name")
        .filter(|part| part.kind() == "type_name")
        .map(|part| unit.name(part))
}

pub(crate) fn procedure_type_keyword(node: Node<'_>) -> Option<&'static str> {
    match node.child(0)?.kind() {
        "double" => Some("real"),
        kind @ ("integer" | "real" | "complex" | "logical" | "character") => Some(kind),
        _ => None,
    }
}

fn use_statement(unit: &Unit, node: Node<'_>) -> Option<Use> {
    let module = child(node, "module_name")?;
    let only = child(node, "included_items");
    let mut cursor = node.walk();
    let nature = node
        .children(&mut cursor)
        .take_while(|item| *item != module)
        .map(|item| item.kind())
        .find(|kind| matches!(*kind, "intrinsic" | "non_intrinsic"));
    let non_intrinsic = nature == Some("non_intrinsic");
    let intrinsic = nature == Some("intrinsic");
    let mut imports = BTreeMap::new();
    let mut imported_at = BTreeMap::new();
    for item in children(only.unwrap_or(node)) {
        let (local, remote) = if item.kind() == "use_alias" {
            let mut names = children(item)
                .into_iter()
                .filter(|part| part.kind() != "comment");
            match (names.next(), names.next()) {
                (Some(local), Some(remote)) => (unit.name(local), unit.name(remote)),
                _ => continue,
            }
        } else if only.is_some() && NAME_NODES.contains(&item.kind()) {
            let name = unit.name(item);
            (name.clone(), name)
        } else {
            continue;
        };
        imported_at.insert(local.clone(), unit.location(item));
        imports.insert(local, remote);
    }
    Some(Use {
        module: unit.name(module),
        location: unit.location(module),
        intrinsic,
        non_intrinsic,
        only: only.is_some(),
        imports,
        imported_at,
    })
}
