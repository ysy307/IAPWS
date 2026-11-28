//! Numerical-kernel constraints and evidence-based real-kind inference.

use std::collections::{BTreeMap, BTreeSet};

use tree_sitter::Node;

use crate::config::{Config, matches};
use crate::diagnostic::Diagnostic;
use crate::semantic::{Project, Symbol, intrinsic_type};
use crate::source::{Unit, child, children};

const PROCEDURES: &[&str] = &["subroutine", "function", "module_procedure"];
const CALLS: &[&str] = &["subroutine_call", "call_expression"];
const ASSIGNMENTS: &[&str] = &[
    "assignment_statement",
    "init_declarator",
    "parameter_assignment",
];
const ARRAY_INTRINSICS: &[&str] = &[
    "cshift",
    "eoshift",
    "matmul",
    "pack",
    "reshape",
    "shape",
    "spread",
    "transpose",
    "unpack",
];
const SAME_KIND_INTRINSICS: &[&str] = &[
    "acos",
    "acosh",
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "bessel_j0",
    "bessel_j1",
    "bessel_y0",
    "bessel_y1",
    "conjg",
    "cos",
    "cosh",
    "dim",
    "epsilon",
    "erf",
    "erfc",
    "erfc_scaled",
    "exp",
    "fraction",
    "gamma",
    "huge",
    "hypot",
    "log",
    "log10",
    "log_gamma",
    "max",
    "maxval",
    "merge",
    "min",
    "minval",
    "mod",
    "modulo",
    "nearest",
    "norm2",
    "product",
    "rrspacing",
    "scale",
    "set_exponent",
    "sign",
    "sin",
    "sinh",
    "spacing",
    "sqrt",
    "sum",
    "tan",
    "tanh",
    "tiny",
];

fn selected(patterns: &[String], value: &str) -> bool {
    patterns.iter().any(|pattern| matches(pattern, value))
}

fn declarator(kind: &str) -> bool {
    matches!(
        kind,
        "init_declarator" | "pointer_init_declarator" | "sized_declarator" | "coarray_declarator"
    )
}

fn type_specification<'a>(symbol: &Symbol<'a>) -> Option<Node<'a>> {
    let mut parent = symbol.node.parent()?;
    while declarator(parent.kind()) {
        parent = parent.parent()?;
    }
    parent.child_by_field_name("type")
}

fn keyword<'a>(unit: &Unit, nodes: &[Node<'a>], names: &[&str]) -> Option<Node<'a>> {
    nodes
        .iter()
        .filter(|node| node.kind() == "keyword_argument")
        .find(|node| {
            node.child_by_field_name("name")
                .is_some_and(|key| names.contains(&unit.name(key).as_str()))
        })
        .and_then(|node| node.child_by_field_name("value"))
}

fn argument<'a>(
    unit: &Unit,
    nodes: &[Node<'a>],
    position: usize,
    names: &[&str],
) -> Option<Node<'a>> {
    keyword(unit, nodes, names).or_else(|| {
        nodes
            .get(position)
            .copied()
            .filter(|node| node.kind() != "keyword_argument")
    })
}

fn arguments(node: Node<'_>) -> Vec<Node<'_>> {
    child(node, "argument_list")
        .map(children)
        .unwrap_or_default()
}

fn external_io(project: &Project<'_>, unit: usize, node: Node<'_>) -> bool {
    let syntax = &project.units[unit];
    match node.kind() {
        "read_statement" | "write_statement" => {
            let target = child(node, "unit_identifier")
                .and_then(|identifier| identifier.named_child(0))
                .or_else(|| keyword(syntax, &children(node), &["unit"]));
            target.and_then(|item| project.expression_category(unit, item)) != Some("character")
        }
        "inquire_statement" => keyword(syntax, &children(node), &["iolength"]).is_none(),
        _ => true,
    }
}

fn procedure_result(project: &Project<'_>, symbol_id: usize) -> Option<usize> {
    let symbol = &project.symbols[symbol_id];
    let scope = symbol
        .procedure
        .or_else(|| project.implementation(symbol_id))?;
    let procedure = &project.scopes[scope];
    let declaration = &project.scopes[procedure.declaration.unwrap_or(scope)];
    let name = procedure.result.as_ref().or(declaration.result.as_ref())?;
    procedure
        .symbols
        .get(name)
        .or_else(|| declaration.symbols.get(name))
        .copied()
}

fn result_symbols(project: &Project<'_>, target_id: usize) -> Vec<usize> {
    let target = &project.symbols[target_id];
    let mut candidates: Vec<_> = if target.kind == "generic" {
        target
            .specifics
            .iter()
            .filter_map(|name| match target.owner {
                Some(owner) => project.member(owner, name),
                None => project.lookup(target.scope, name),
            })
            .collect()
    } else {
        vec![target_id]
    };
    if target.kind == "generic" && candidates.len() != target.specifics.len() {
        return Vec::new();
    }
    candidates.retain_mut(|candidate| {
        let symbol = &project.symbols[*candidate];
        if symbol.kind == "binding" {
            let Some(id) = project.binding_target(*candidate) else {
                return false;
            };
            *candidate = id;
        }
        let symbol = &project.symbols[*candidate];
        if let Some(interface) = symbol.interface_name.as_deref()
            && let Some(id) = project.lookup_kind(
                symbol.scope,
                interface,
                &["procedure", "abstract_interface"],
            )
        {
            *candidate = id;
        }
        if let Some(result) = procedure_result(project, *candidate) {
            *candidate = result;
            true
        } else {
            false
        }
    });
    candidates
}

fn intrinsic_name<'a>(project: &'a Project<'_>, unit: usize, node: Node<'_>) -> Option<&'a str> {
    let id = node
        .named_child(0)
        .and_then(|target| project.expression(unit, target))
        .or_else(|| project.expression(unit, node))?;
    let symbol = &project.symbols[id];
    (symbol.kind == "intrinsic" || symbol.attrs.contains("intrinsic"))
        .then_some(symbol.name.as_str())
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Category {
    Real,
    Complex,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct RealType {
    category: Category,
    width: Option<usize>,
}

struct Precision<'a, 'p> {
    project: &'p Project<'a>,
    widths: &'p BTreeMap<String, usize>,
    expressions: BTreeMap<(usize, usize), Option<RealType>>,
    parameters: BTreeMap<usize, (usize, Node<'a>)>,
    active_parameters: BTreeSet<usize>,
}

impl<'a, 'p> Precision<'a, 'p> {
    fn new(project: &'p Project<'a>, widths: &'p BTreeMap<String, usize>) -> Self {
        let mut parameters = BTreeMap::new();
        for scope in &project.scopes {
            let unit = &project.units[scope.unit];
            for statement in children(scope.node)
                .into_iter()
                .filter(|node| node.kind() == "parameter_statement")
            {
                for assignment in children(statement) {
                    if let (Some(name), Some(value)) =
                        (assignment.named_child(0), assignment.named_child(1))
                        && let Some(id) = scope.symbols.get(&unit.name(name))
                    {
                        parameters.insert(*id, (scope.unit, value));
                    }
                }
            }
        }
        Self {
            project,
            widths,
            expressions: BTreeMap::new(),
            parameters,
            active_parameters: BTreeSet::new(),
        }
    }

    fn width(&mut self, unit: usize, value: Option<Node<'a>>) -> Option<usize> {
        let node = value?;
        match node.kind() {
            "parenthesized_expression" => self.width(unit, node.named_child(0)),
            "call_expression" if intrinsic_name(self.project, unit, node) == Some("kind") => {
                let first = argument(&self.project.units[unit], &arguments(node), 0, &["x"]);
                self.expression_type(unit, first)?.width
            }
            "identifier" => {
                let id = self.project.expression(unit, node)?;
                let symbol = &self.project.symbols[id];
                if let Some(origin) = &symbol.external {
                    return self
                        .widths
                        .get(&format!("{origin}::{}", symbol.name))
                        .copied();
                }
                if symbol.kind != "constant" || !self.active_parameters.insert(id) {
                    return None;
                }
                let initializer = symbol
                    .node
                    .parent()
                    .filter(|parent| parent.kind() == "init_declarator")
                    .and_then(|parent| parent.child_by_field_name("right"))
                    .map(|value| (self.project.scopes[symbol.scope].unit, value))
                    .or_else(|| self.parameters.get(&id).copied());
                let size = initializer.and_then(|(origin_unit, expression)| {
                    self.width(origin_unit, Some(expression))
                });
                self.active_parameters.remove(&id);
                size
            }
            _ => None,
        }
    }

    fn declared_type(&mut self, unit: usize, specification: Option<Node<'a>>) -> Option<RealType> {
        let specification = specification?;
        let category = match intrinsic_type(Some(specification))? {
            "real" => Category::Real,
            "complex" => Category::Complex,
            _ => return None,
        };
        let parameter = specification
            .child_by_field_name("kind")
            .and_then(|selector| {
                argument(&self.project.units[unit], &children(selector), 0, &["kind"])
            });
        Some(RealType {
            category,
            width: self.width(unit, parameter),
        })
    }

    fn symbol_type(&mut self, id: usize) -> Option<RealType> {
        let symbol = &self.project.symbols[id];
        let unit = self.project.scopes[symbol.scope].unit;
        self.declared_type(unit, type_specification(symbol))
            .or_else(|| {
                let category = match symbol.intrinsic.as_deref()? {
                    "real" => Category::Real,
                    "complex" => Category::Complex,
                    _ => return None,
                };
                Some(RealType {
                    category,
                    width: None,
                })
            })
    }

    fn expression_type(&mut self, unit: usize, value: Option<Node<'a>>) -> Option<RealType> {
        let node = value?;
        let key = (unit, node.id());
        if let Some(result) = self.expressions.get(&key) {
            return *result;
        }
        // Named kind constants may themselves query the type of an expression.
        self.expressions.insert(key, None);
        let result = self.infer(unit, node);
        self.expressions.insert(key, result);
        result
    }

    fn integer(&self, unit: usize, node: Node<'a>) -> bool {
        match node.kind() {
            "math_expression" => children(node)
                .into_iter()
                .all(|item| self.integer(unit, item)),
            "parenthesized_expression" | "unary_expression" => node
                .named_child(0)
                .is_some_and(|item| self.integer(unit, item)),
            _ => self.project.expression_category(unit, node) == Some("integer"),
        }
    }

    fn arithmetic(&mut self, unit: usize, left: Node<'a>, right: Node<'a>) -> Option<RealType> {
        let lhs = self.expression_type(unit, Some(left));
        let rhs = self.expression_type(unit, Some(right));
        match (lhs, rhs) {
            (Some(a), Some(b)) => Some(RealType {
                category: a.category.max(b.category),
                width: a.width.zip(b.width).map(|(a, b)| a.max(b)),
            }),
            (Some(a), None) if self.integer(unit, right) => Some(a),
            (None, Some(b)) if self.integer(unit, left) => Some(b),
            _ => None,
        }
    }

    fn call_type(&mut self, unit: usize, node: Node<'a>) -> Option<RealType> {
        let project = self.project;
        let Some(name) = intrinsic_name(project, unit, node) else {
            let id = project.expression(unit, node)?;
            let results = result_symbols(project, id);
            return if results.is_empty() {
                self.symbol_type(id)
            } else {
                self.uniform_symbols(&results)
            };
        };
        let args = arguments(node);
        let syntax = &project.units[unit];
        let first = argument(
            syntax,
            &args,
            0,
            &[
                "a", "a1", "array", "array_a", "matrix", "matrix_a", "source", "t", "tsource", "x",
                "vector_a", "z",
            ],
        );
        let first_type = self.expression_type(unit, first);
        match name {
            "real" | "cmplx" | "aint" | "anint" => {
                let parameter = argument(
                    syntax,
                    &args,
                    if name == "cmplx" { 2 } else { 1 },
                    &["kind"],
                );
                let category = if name == "cmplx" {
                    Category::Complex
                } else {
                    Category::Real
                };
                let width = if parameter.is_none()
                    && (matches!(name, "aint" | "anint")
                        || (name == "real"
                            && first_type.is_some_and(|kind| kind.category == Category::Complex)))
                {
                    first_type.and_then(|kind| kind.width)
                } else {
                    self.width(unit, parameter)
                };
                Some(RealType { category, width })
            }
            "abs" | "aimag" => Some(RealType {
                category: Category::Real,
                width: first_type?.width,
            }),
            "dble" => Some(RealType {
                category: Category::Real,
                width: None,
            }),
            "dot_product" | "matmul" => {
                let rhs = argument(syntax, &args, 1, &["array_b", "matrix_b", "vector_b"])?;
                self.arithmetic(unit, first?, rhs)
            }
            "shape" => None,
            _ if SAME_KIND_INTRINSICS.contains(&name) || ARRAY_INTRINSICS.contains(&name) => {
                first_type
            }
            _ => None,
        }
    }

    fn uniform_symbols(&mut self, ids: &[usize]) -> Option<RealType> {
        let mut kinds = ids.iter().map(|id| self.symbol_type(*id));
        let first = kinds.next()??;
        kinds.all(|value| value == Some(first)).then_some(first)
    }

    fn infer(&mut self, unit: usize, node: Node<'a>) -> Option<RealType> {
        match node.kind() {
            "number_literal" if self.project.expression_category(unit, node) == Some("real") => {
                Some(RealType {
                    category: Category::Real,
                    width: self.width(unit, node.child_by_field_name("kind")),
                })
            }
            "math_expression" => self.arithmetic(
                unit,
                node.child_by_field_name("left")?,
                node.child_by_field_name("right")?,
            ),
            "parenthesized_expression" | "unary_expression" => {
                self.expression_type(unit, node.named_child(0))
            }
            "complex_literal" => {
                let width = self
                    .arithmetic(unit, node.named_child(0)?, node.named_child(1)?)
                    .and_then(|kind| kind.width);
                Some(RealType {
                    category: Category::Complex,
                    width,
                })
            }
            "array_literal" | "implied_do_loop_expression" => {
                if let Some(specification) = child(node, "intrinsic_type") {
                    return self.declared_type(unit, Some(specification));
                }
                let mut kinds = children(node)
                    .into_iter()
                    .filter(|item| item.kind() != "loop_control_expression")
                    .map(|item| self.expression_type(unit, Some(item)));
                let first = kinds.next()??;
                kinds.all(|kind| kind == Some(first)).then_some(first)
            }
            "call_expression" => self.call_type(unit, node),
            _ => self.symbol_type(self.project.expression(unit, node)?),
        }
    }

    fn mixed(&mut self, unit: usize, node: Node<'a>) -> BTreeSet<usize> {
        let operands = match node.kind() {
            "math_expression" | "relational_expression" => vec![
                node.child_by_field_name("left"),
                node.child_by_field_name("right"),
            ],
            "complex_literal" => children(node).into_iter().map(Some).collect(),
            "array_literal" if child(node, "intrinsic_type").is_none() => {
                children(node).into_iter().map(Some).collect()
            }
            "call_expression"
                if intrinsic_name(self.project, unit, node).is_some_and(|name| {
                    ["min", "max", "matmul", "dot_product", "merge"].contains(&name)
                }) =>
            {
                arguments(node)
                    .into_iter()
                    .map(|item| {
                        if item.kind() == "keyword_argument" {
                            item.child_by_field_name("value")
                        } else {
                            Some(item)
                        }
                    })
                    .collect()
            }
            _ => return BTreeSet::new(),
        };
        operands
            .into_iter()
            .filter_map(|operand| self.expression_type(unit, operand)?.width)
            .collect()
    }

    fn narrowing(&mut self, unit: usize, node: Node<'a>) -> Option<String> {
        let (target, value) = match node.kind() {
            "init_declarator" => {
                let kind = node
                    .parent()
                    .and_then(|parent| parent.child_by_field_name("type"));
                (
                    self.declared_type(unit, kind),
                    node.child_by_field_name("right"),
                )
            }
            "parameter_assignment" => (
                self.expression_type(unit, node.named_child(0)),
                node.named_child(1),
            ),
            _ => (
                self.expression_type(unit, node.child_by_field_name("left")),
                node.child_by_field_name("right"),
            ),
        };
        let target = target?;
        let value = self.expression_type(unit, value)?;
        if target.category == Category::Real && value.category == Category::Complex {
            return Some("Implicit assignment discards the imaginary component".into());
        }
        let (source, destination) = (value.width?, target.width?);
        (source > destination).then(|| {
            format!("Implicit assignment narrows floating-point precision from {source} to {destination} bits")
        })
    }
}

fn polymorphic(project: &Project<'_>, unit: usize, node: Node<'_>) -> bool {
    let Some(id) = project.expression(unit, node) else {
        return false;
    };
    let symbol = &project.symbols[id];
    if symbol.kind == "type" && node.kind() == "identifier" {
        let ancestors = std::iter::successors(node.parent(), |ancestor| ancestor.parent());
        for ancestor in ancestors.filter(|ancestor| ancestor.kind() == "type_statement") {
            if ancestor
                .child_by_field_name("type")
                .and_then(|guard| project.expression(unit, guard))
                == Some(id)
            {
                return ancestor
                    .child(0)
                    .is_some_and(|token| token.kind() == "class");
            }
        }
    }
    let targets = result_symbols(project, id);
    let symbols = if targets.is_empty() {
        vec![id]
    } else {
        targets
    };
    symbols.iter().any(|id| {
        type_specification(&project.symbols[*id]).is_some_and(|spec| {
            spec.kind() == "derived_type"
                && spec.child(0).is_some_and(|token| token.kind() == "class")
        })
    })
}

fn dynamic_dispatch(project: &Project<'_>, unit: usize, node: Node<'_>) -> bool {
    let Some(target) = node.named_child(0) else {
        return false;
    };
    let Some(id) = project
        .expression(unit, node)
        .or_else(|| project.expression(unit, target))
    else {
        return false;
    };
    let symbol = &project.symbols[id];
    if symbol.kind == "procedure" && symbol.attrs.contains("pointer") {
        return true;
    }
    let Some(owner) = symbol.owner else {
        return false;
    };
    if !matches!(symbol.kind, "binding" | "generic")
        || symbol.attrs.contains("non_overridable")
        || target.kind() != "derived_type_member_expression"
        || !target
            .named_child(0)
            .is_some_and(|item| polymorphic(project, unit, item))
    {
        return false;
    }
    symbol.kind == "binding"
        || symbol.specifics.iter().any(|name| {
            project
                .member(owner, name)
                .is_some_and(|id| !project.symbols[id].attrs.contains("non_overridable"))
        })
}

fn array_symbol(symbol: &Symbol<'_>) -> bool {
    symbol
        .attrs
        .iter()
        .any(|attr| attr == "dimension" || attr.starts_with("dimension("))
        || std::iter::successors(symbol.node.parent(), |node| node.parent())
            .take_while(|node| declarator(node.kind()))
            .any(|node| node.kind() == "sized_declarator")
}

fn array_expression(project: &Project<'_>, unit: usize, node: Node<'_>) -> bool {
    match node.kind() {
        "array_literal" => true,
        "math_expression" | "parenthesized_expression" | "unary_expression" => children(node)
            .into_iter()
            .any(|item| array_expression(project, unit, item)),
        "call_expression" => {
            let Some(id) = project.expression(unit, node) else {
                return false;
            };
            if array_symbol(&project.symbols[id]) {
                return arguments(node).into_iter().any(|item| {
                    item.kind() == "extent_specifier" || array_expression(project, unit, item)
                });
            }
            result_symbols(project, id)
                .into_iter()
                .any(|id| array_symbol(&project.symbols[id]))
                || intrinsic_name(project, unit, node)
                    .is_some_and(|name| ARRAY_INTRINSICS.contains(&name))
        }
        _ => project
            .expression(unit, node)
            .is_some_and(|id| array_symbol(&project.symbols[id])),
    }
}

fn preserves_shape(
    project: &Project<'_>,
    unit: usize,
    target: Node<'_>,
    expression: Node<'_>,
) -> bool {
    let syntax = &project.units[unit];
    let same_symbol = project.expression(unit, target) == project.expression(unit, expression);
    if same_symbol && syntax.name(target) == syntax.name(expression) {
        return true;
    }
    matches!(
        expression.kind(),
        "math_expression" | "parenthesized_expression" | "unary_expression"
    ) && children(expression)
        .into_iter()
        .any(|item| preserves_shape(project, unit, target, item))
}

fn allocation(project: &Project<'_>, unit: usize, node: Node<'_>) -> Option<&'static str> {
    match node.kind() {
        "allocate_statement" => Some("ALLOCATE is executed inside a designated hot loop"),
        "call_expression" => {
            let target = project.expression(unit, node)?;
            result_symbols(project, target)
                .into_iter()
                .any(|id| project.symbols[id].attrs.contains("allocatable"))
                .then_some("Allocatable function result may allocate inside a designated hot loop")
        }
        "assignment_statement" => {
            let lhs = node.child_by_field_name("left")?;
            let rhs = node.child_by_field_name("right")?;
            if !matches!(lhs.kind(), "identifier" | "derived_type_member_expression") {
                return None;
            }
            let target = &project.symbols[project.expression(unit, lhs)?];
            if !target.attrs.contains("allocatable") {
                return None;
            }
            if array_symbol(target)
                && array_expression(project, unit, rhs)
                && !preserves_shape(project, unit, lhs, rhs)
            {
                return Some(
                    "Whole allocatable array assignment may reallocate inside a designated hot loop",
                );
            }
            let origin = &project.units[project.scopes[target.scope].unit];
            let selector =
                type_specification(target).and_then(|spec| spec.child_by_field_name("kind"));
            let size =
                selector.and_then(|parameter| argument(origin, &children(parameter), 0, &["len"]));
            (target.intrinsic.as_deref() == Some("character")
                && size.is_some_and(|node| origin.name(node) == ":")
                && !preserves_shape(project, unit, lhs, rhs))
            .then_some(
                "Deferred-length character assignment may reallocate inside a designated hot loop",
            )
        }
        _ => None,
    }
}

pub fn check(project: &Project<'_>, config: &Config) -> Vec<Diagnostic> {
    let policy = &config.numerical;
    let mut precision = Precision::new(project, &policy.kinds);
    let mut diagnostics = Vec::new();
    for (id, scope) in project.scopes.iter().enumerate() {
        // BLOCK declarations and bodies are visited with their containing procedure.
        if scope.kind == "block" {
            continue;
        }
        let declaration = &project.scopes[scope.declaration.unwrap_or(id)];
        let names = [&scope.qualified, &declaration.qualified];
        let selects = |patterns: &[String]| names.iter().any(|name| selected(patterns, name));
        let kernel =
            PROCEDURES.contains(&scope.kind) && !scope.interface_body && selects(&policy.kernels);
        let hot = selects(&policy.hot_loops);
        if kernel
            && (declaration.attrs.contains("impure")
                || (!declaration.attrs.contains("pure")
                    && !declaration.attrs.contains("elemental")))
        {
            let syntax = &project.units[declaration.unit];
            let diagnostic = Diagnostic::new(
                "NUM001",
                "Numerical kernel lacks PURE or implicitly pure ELEMENTAL",
                syntax.location(declaration.header),
            )
            .with_subject(&declaration.qualified);
            diagnostics.push(if selects(&policy.purity_warnings) {
                diagnostic.warning()
            } else {
                diagnostic
            });
        }

        let unit = scope.unit;
        let syntax = &project.units[unit];
        let mut pending = vec![(scope.node, false)];
        while let Some((node, in_body)) = pending.pop() {
            let loop_node = hot && selected(&policy.hot_loop_nodes, node.kind());
            let repeated = in_body || (loop_node && node.kind() != "do_loop");
            let mut emit = |code, message: String, warning| {
                let item = Diagnostic::new(code, message, syntax.location(node))
                    .with_subject(&declaration.qualified);
                diagnostics.push(if warning { item.warning() } else { item });
            };
            if kernel
                && selected(&policy.io_statements, node.kind())
                && external_io(project, unit, node)
            {
                emit("NUM002", "External I/O in a numerical kernel".into(), false);
            }
            let widths = precision.mixed(unit, node);
            if widths.len() > 1 {
                let sizes = widths
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(", ");
                emit(
                    "NUM003",
                    format!("Mixed floating-point kinds with storage widths {sizes} bits"),
                    true,
                );
            }
            if ASSIGNMENTS.contains(&node.kind())
                && let Some(message) = precision.narrowing(unit, node)
            {
                emit("NUM004", message, false);
            }
            if repeated {
                if CALLS.contains(&node.kind()) && dynamic_dispatch(project, unit, node) {
                    emit(
                        "NUM005",
                        "Dynamic procedure dispatch inside a designated hot loop".into(),
                        true,
                    );
                }
                if let Some(message) = allocation(project, unit, node) {
                    emit("NUM006", message.into(), true);
                }
            }
            pending.extend(
                children(node)
                    .into_iter()
                    .rev()
                    .filter(|item| {
                        !PROCEDURES.contains(&item.kind())
                            && !matches!(item.kind(), "internal_procedures" | "interface")
                    })
                    .map(|item| {
                        (
                            item,
                            repeated || (loop_node && item.kind() != "do_statement"),
                        )
                    }),
            );
        }
    }
    diagnostics
}
