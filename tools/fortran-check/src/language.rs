//! Check source conventions and encapsulation before applying configured exceptions.

use crate::config::Config;
use crate::diagnostic::Diagnostic;
use crate::semantic::{Project, attributes, declaration_name, procedure_type_keyword, spelling};
use crate::source::{Unit, child, children, fields, walk};
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::LazyLock;
use tree_sitter::Node;

static REAL_LITERAL: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?i)^(?:(?:\d+\.\d*|\.\d+)(?:[edq][+-]?\d+)?|\d+[edq][+-]?\d+)$")
        .expect("valid Fortran real literal regular expression")
});

const LITERALS: &[&str] = &["number_literal", "boolean_literal", "string_literal"];

pub fn check(project: &Project<'_>, config: &Config) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    let binding_targets = exposed_implementations(project);
    for (id, scope) in project.scopes.iter().enumerate() {
        check_scope_language(project, id, config, &mut diagnostics);
        if matches!(scope.kind, "module" | "submodule") {
            for &id in scope.symbols.values() {
                let symbol = &project.symbols[id];
                if symbol.kind == "variable"
                    || symbol.kind == "procedure" && symbol.attrs.contains("pointer")
                {
                    diagnostics.push(symbol_issue(
                        project,
                        "OOP003",
                        "mutable module variable is prohibited",
                        id,
                    ));
                }
            }
        }
        if scope.kind == "module" {
            check_public_api(project, id, &binding_targets, &mut diagnostics);
        } else if matches!(scope.kind, "subroutine" | "function") && scope.parent.is_none() {
            diagnostics.push(
                Diagnostic::new(
                    "OOP005",
                    "external procedure requires a documented exception",
                    project.unit(id).location(scope.header),
                )
                .with_subject(&scope.qualified),
            );
        }
    }
    let mut type_definitions = BTreeMap::new();
    for &owner in &project.types {
        check_type(project, owner, &mut diagnostics);
        let owner = &project.symbols[owner];
        if let Some(definition) = owner.node.parent().and_then(|node| node.parent()) {
            type_definitions.insert(
                (project.scopes[owner.scope].unit, definition.id()),
                &owner.qualified,
            );
        }
    }
    let banned: BTreeSet<_> = config
        .language
        .banned_statements
        .iter()
        .map(|name| statement_key(name))
        .collect();
    for (index, unit) in project.units.iter().enumerate() {
        for node in walk(unit.tree.root_node()) {
            if let Some(location) = bad_kind(project, index, node, config) {
                diagnostics.push(
                    Diagnostic::new(
                        "STD004",
                        "use a named kind constant",
                        unit.location(location),
                    )
                    .with_subject(context(
                        project,
                        &type_definitions,
                        index,
                        location,
                    )),
                );
            }
            if let Some(name) = banned_statement(unit, node, &banned) {
                diagnostics.push(
                    Diagnostic::new(
                        "STD006",
                        format!("{name} is prohibited"),
                        unit.location(node),
                    )
                    .with_subject(context(
                        project,
                        &type_definitions,
                        index,
                        node,
                    )),
                );
            }
            if node.kind() == "derived_type" && child(node, "unlimited_polymorphic").is_some() {
                check_polymorphism(project, &type_definitions, index, node, &mut diagnostics);
            }
        }
    }
    diagnostics
}

fn check_scope_language(
    project: &Project<'_>,
    id: usize,
    config: &Config,
    diagnostics: &mut Vec<Diagnostic>,
) {
    let scope = &project.scopes[id];
    let declaration = &project.scopes[scope.declaration.unwrap_or(id)];
    let unit = project.unit(id);
    let statements = children(scope.node);
    let mut emit = |code, message: String, node, name: &str| {
        let subject = entity_subject(&declaration.qualified, name);
        diagnostics.push(Diagnostic::new(code, message, unit.location(node)).with_subject(subject));
    };
    let mut implicit = statements
        .iter()
        .copied()
        .filter(|node| node.kind() == "implicit_statement");
    let first = implicit.next();
    let inherited =
        scope.kind == "module_procedure" && scope.declaration.is_some() && first.is_none();
    let require_external = config.language.require_implicit_external;
    if config
        .language
        .implicit_scopes
        .iter()
        .any(|candidate| candidate.eq_ignore_ascii_case(scope.kind))
        && !inherited
        && !(implicit.next().is_none()
            && first.is_some_and(|node| implicit_none(node, require_external)))
    {
        let declaration = if require_external {
            "IMPLICIT NONE (TYPE, EXTERNAL)"
        } else {
            "IMPLICIT NONE"
        };
        emit(
            "STD001",
            format!("require {declaration} in this scoping unit"),
            first.unwrap_or(scope.header),
            "",
        );
    }
    for node in statements
        .iter()
        .copied()
        .filter(|node| node.kind() == "use_statement")
    {
        if child(node, "included_items").is_none()
            && let Some(module) = child(node, "module_name")
        {
            emit(
                "C121",
                format!("USE of '{}' requires ONLY", unit.name(module)),
                node,
                "",
            );
        }
    }

    let mut saved = BTreeSet::new();
    let mut save_all = false;
    let mut modifications: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for &node in &statements {
        match node.kind() {
            "save_statement" => {
                let names = significant_children(node);
                save_all |= names.is_empty();
                saved.extend(
                    names
                        .into_iter()
                        .filter(|name| name.kind() == "identifier")
                        .map(|name| unit.name(name)),
                );
            }
            "parameter_statement" => {
                for item in children(node)
                    .into_iter()
                    .filter(|node| node.kind() == "parameter_assignment")
                {
                    if let Some(name) = entity(item) {
                        modifications
                            .entry(unit.name(name))
                            .or_default()
                            .insert("parameter".into());
                    }
                }
            }
            "variable_modification" => {
                let attrs = attributes(unit, node);
                for item in fields(node, "declarator").into_iter().filter_map(entity) {
                    modifications
                        .entry(unit.name(item))
                        .or_default()
                        .extend(attrs.iter().cloned());
                }
            }
            _ => {}
        }
    }
    let mut reported_save = BTreeSet::new();
    for &node in &statements {
        if !matches!(
            node.kind(),
            "variable_declaration" | "variable_modification"
        ) {
            continue;
        }
        let attrs = attributes(unit, node);
        let procedure = node
            .child_by_field_name("type")
            .filter(|node| node.kind() == "procedure");
        let implicit_interface = procedure.is_some_and(|node| {
            child(node, "procedure_interface").is_none_or(|interface| {
                procedure_type_keyword(interface).is_some()
                    && project.expression(scope.unit, interface).is_none()
            })
        });
        for item in fields(node, "declarator") {
            let Some(name_node) = entity(item) else {
                continue;
            };
            let name = unit.name(name_node);
            if attrs.contains("intent(out)") {
                emit(
                    "STD002",
                    format!("INTENT(OUT) requires a documented exception for '{name}'"),
                    name_node,
                    &name,
                );
            }
            if attrs.contains("external") || implicit_interface {
                emit(
                    "STD003",
                    format!("'{name}' requires an explicit procedure interface"),
                    name_node,
                    &name,
                );
            }
            if node.kind() == "variable_modification" {
                continue;
            }
            let has = |attribute| {
                attrs.contains(attribute)
                    || modifications
                        .get(&name)
                        .is_some_and(|attrs| attrs.contains(attribute))
            };
            let explicitly_saved = save_all || saved.contains(&name) || has("save");
            if explicitly_saved {
                saved.insert(name.clone());
            }
            if !explicitly_saved
                && !has("parameter")
                && !has("external")
                && !has("intrinsic")
                && (procedure.is_none() || has("pointer"))
                && (matches!(scope.kind, "module" | "submodule" | "program")
                    || matches!(item.kind(), "init_declarator" | "pointer_init_declarator"))
            {
                emit(
                    "STD005",
                    format!("'{name}' has implicit SAVE"),
                    name_node,
                    &name,
                );
                reported_save.insert(name);
            }
        }
    }

    if !save_all {
        for data in statements
            .into_iter()
            .filter(|node| node.kind() == "data_statement")
        {
            for data_set in children(data)
                .into_iter()
                .filter(|node| node.kind() == "data_set")
            {
                let mut targets = Vec::new();
                data_targets(data_set, &mut targets);
                for name_node in targets {
                    let name = unit.name(name_node);
                    if !saved.contains(&name) && reported_save.insert(name.clone()) {
                        emit(
                            "STD005",
                            format!("DATA initialization gives '{name}' implicit SAVE"),
                            name_node,
                            &name,
                        );
                    }
                }
            }
        }
    }
}

fn implicit_none(node: Node<'_>, require_external: bool) -> bool {
    let mut cursor = node.walk();
    let parts: BTreeSet<_> = node.children(&mut cursor).map(|part| part.kind()).collect();
    let external = parts.contains("external");
    parts.contains("none") && (parts.contains("type") || !external) && (!require_external || external)
}

fn entity(node: Node<'_>) -> Option<Node<'_>> {
    let node = declaration_name(node);
    matches!(node.kind(), "identifier" | "name" | "method_name").then_some(node)
}

fn data_targets<'a>(node: Node<'a>, targets: &mut Vec<Node<'a>>) {
    if matches!(node.kind(), "data_set" | "implied_do_loop_expression") {
        for part in significant_children(node) {
            if !matches!(part.kind(), "data_value" | "loop_control_expression") {
                data_targets(part, targets);
            }
        }
    } else if let Some(name) = entity(node) {
        targets.push(name);
    }
}

fn significant_children(node: Node<'_>) -> Vec<Node<'_>> {
    children(node)
        .into_iter()
        .filter(|part| part.kind() != "comment")
        .collect()
}

fn kind_argument<'a>(unit: &Unit, arguments: &[Node<'a>], position: usize) -> Option<Node<'a>> {
    for &argument in arguments {
        if argument.kind() == "keyword_argument"
            && argument
                .child_by_field_name("name")
                .is_some_and(|name| unit.name(name) == "kind")
        {
            return argument.child_by_field_name("value");
        }
    }
    arguments
        .get(position)
        .copied()
        .filter(|argument| argument.kind() != "keyword_argument")
}

fn conversion_kind<'a>(project: &Project<'_>, unit: usize, node: Node<'a>) -> Option<Node<'a>> {
    let children = significant_children(node);
    let target = children
        .first()
        .filter(|node| node.kind() == "identifier")?;
    let args = children
        .get(1)
        .filter(|node| node.kind() == "argument_list")?;
    let source = &project.units[unit];
    let name = source.name(*target);
    let mut position = match name.as_str() {
        "achar" | "aint" | "anint" | "ceiling" | "char" | "floor" | "iachar" | "ichar" | "int"
        | "len" | "len_trim" | "logical" | "maskl" | "maskr" | "nint" | "real" | "shape"
        | "storage_size" => 1,
        "cmplx" | "count" | "lbound" | "size" | "ubound" => 2,
        "index" | "maxloc" | "minloc" | "scan" | "verify" => 3,
        "findloc" => 4,
        _ => return None,
    };
    if let Some(symbol) = project.expression(unit, *target) {
        let symbol = &project.symbols[symbol];
        if symbol.kind != "intrinsic" && !symbol.attrs.contains("intrinsic") {
            return None;
        }
    }
    let arguments = significant_children(*args);
    let first_optional = match name.as_str() {
        "findloc" => arguments.get(2),
        "maxloc" | "minloc" => arguments.get(1),
        _ => None,
    };
    if first_optional
        .is_some_and(|part| project.expression_category(unit, *part) == Some("logical"))
    {
        position -= 1;
    }
    kind_argument(source, &arguments, position)
}

fn bad_kind<'a>(
    project: &Project<'_>,
    unit: usize,
    node: Node<'a>,
    config: &Config,
) -> Option<Node<'a>> {
    let require_named = |category: &str| {
        config
            .language
            .named_kinds
            .iter()
            .any(|name| name.eq_ignore_ascii_case(category))
    };
    let (required, selector) = match node.kind() {
        "intrinsic_type" => {
            let category = node.child(0)?.kind();
            let kind = node.child_by_field_name("kind").and_then(|part| {
                kind_argument(
                    &project.units[unit],
                    &significant_children(part),
                    usize::from(category == "character"),
                )
            });
            (require_named(category), kind)
        }
        "call_expression" => (false, conversion_kind(project, unit, node)),
        kind if LITERALS.contains(&kind) => {
            let selector = node.child_by_field_name("kind");
            let real = kind == "number_literal"
                && selector.is_none()
                && require_named("real")
                && REAL_LITERAL.is_match(project.units[unit].text(node));
            (real, selector)
        }
        _ => return None,
    };
    match selector {
        Some(selector) if selector.kind() != "identifier" => Some(selector),
        None if required => Some(node),
        _ => None,
    }
}

fn statement_key(name: &str) -> String {
    name.chars()
        .filter(char::is_ascii_alphanumeric)
        .map(|character| character.to_ascii_lowercase())
        .collect()
}

fn banned_statement(unit: &Unit, node: Node<'_>, banned: &BTreeSet<String>) -> Option<String> {
    if banned.is_empty() {
        return None;
    }
    let kind = node.kind();
    let category = if kind == "keyword_statement" {
        let token = node.child(0)?;
        let first = token.kind();
        if matches!(first, "go" | "goto") {
            "go to"
        } else {
            first
        }
    } else {
        kind.strip_suffix("_statement")?
    };
    if banned.contains(&statement_key(category)) || banned.contains(&statement_key(kind)) {
        let name = if category == "format" {
            spelling(unit, node.child(0)?)
        } else {
            category.replace('_', " ")
        };
        Some(name.to_ascii_uppercase())
    } else {
        None
    }
}

fn exposed_implementations(project: &Project<'_>) -> BTreeSet<usize> {
    project
        .types
        .iter()
        .flat_map(|&id| project.symbols[id].members.values().copied())
        .filter(|&id| {
            let member = &project.symbols[id];
            matches!(member.kind, "binding" | "final") && !member.attrs.contains("deferred")
        })
        .filter_map(|id| project.binding_target(id))
        .collect()
}

fn symbol_issue(
    project: &Project<'_>,
    code: &str,
    message: impl Into<String>,
    id: usize,
) -> Diagnostic {
    let symbol = &project.symbols[id];
    Diagnostic::new(
        code,
        message,
        project.unit(symbol.scope).location(symbol.node),
    )
    .with_subject(&symbol.qualified)
}

fn check_type(project: &Project<'_>, owner: usize, diagnostics: &mut Vec<Diagnostic>) {
    let definition = &project.symbols[owner];
    for &id in definition.members.values() {
        let member = &project.symbols[id];
        if matches!(member.kind, "variable" | "procedure")
            && !member.attrs.contains("kind")
            && !member.attrs.contains("len")
        {
            // Component-access statements are permitted for module types.
            if project.scopes[definition.scope].kind == "module" && member.is_public {
                diagnostics.push(symbol_issue(
                    project,
                    "OOP002",
                    "derived-type component must be private",
                    id,
                ));
            }
            if member.kind == "variable" && member.attrs.contains("pointer") {
                diagnostics.push(symbol_issue(
                    project,
                    "OOP010",
                    "data pointer component requires an exception",
                    id,
                ));
            }
            if member.kind == "procedure" && member.interface_name.is_none() {
                diagnostics.push(symbol_issue(
                    project,
                    "STD003",
                    format!("'{}' requires an explicit procedure interface", member.name),
                    id,
                ));
            }
        } else if matches!(member.kind, "binding" | "generic") {
            if !member.attrs.contains("public") && !member.attrs.contains("private") {
                diagnostics.push(symbol_issue(
                    project,
                    "OOP004",
                    "specify PUBLIC or PRIVATE on this binding",
                    id,
                ));
            }
            if member.kind == "generic" && member.is_public {
                check_generic(project, id, diagnostics);
            }
        }
        if member.attrs.contains("nopass") {
            diagnostics.push(symbol_issue(
                project,
                "OOP009",
                "NOPASS requires a documented exception",
                id,
            ));
        }
    }
}

fn check_generic(project: &Project<'_>, id: usize, diagnostics: &mut Vec<Diagnostic>) {
    let generic = &project.symbols[id];
    if generic.specifics.is_empty() {
        diagnostics.push(symbol_issue(
            project,
            "OOP011",
            "public generic has no specific procedures",
            id,
        ));
        return;
    }
    let mut seen = BTreeSet::new();
    for name in &generic.specifics {
        let valid = if let Some(owner) = generic.owner {
            project.member(owner, name).is_some_and(|target| {
                let target = &project.symbols[target];
                target.kind == "binding" && (target.is_public || target.scope == generic.scope)
            })
        } else {
            project
                .lookup_kind(generic.scope, name, &["procedure"])
                .is_some_and(|target| {
                    let target = &project.symbols[target];
                    matches!(target.kind, "procedure" | "external")
                        && !target.attrs.contains("external")
                        && !target.attrs.contains("pointer")
                })
        };
        if !seen.insert(name) || !valid {
            diagnostics.push(symbol_issue(
                project,
                "OOP011",
                format!(
                    "public generic '{}' has invalid specific '{name}'",
                    generic.name
                ),
                id,
            ));
        }
    }
}

fn check_public_api(
    project: &Project<'_>,
    id: usize,
    bound: &BTreeSet<usize>,
    diagnostics: &mut Vec<Diagnostic>,
) {
    let scope = &project.scopes[id];
    let unit = project.unit(id);
    if !scope.default_private {
        diagnostics.push(
            Diagnostic::new(
                "OOP001",
                "module default accessibility must be PRIVATE",
                unit.location(scope.header),
            )
            .with_subject(&scope.qualified),
        );
    }
    let mut exports: Vec<_> = scope
        .symbols
        .iter()
        .chain(&scope.generics)
        .filter(|&(_, &id)| {
            let symbol = &project.symbols[id];
            symbol.is_public && matches!(symbol.kind, "procedure" | "generic")
        })
        .map(|(name, &id)| (name, id))
        .collect();
    for (name, &public) in &scope.explicit_access {
        if public
            && !scope.symbols.contains_key(name)
            && !scope.generics.contains_key(name)
            && let Some(symbol) = project.lookup(id, name)
        {
            exports.push((name, symbol));
        }
    }
    for (name, symbol_id) in exports {
        let symbol = &project.symbols[symbol_id];
        let location = scope.explicit_access_nodes.get(name).map_or_else(
            || project.unit(symbol.scope).location(symbol.node),
            |&node| unit.location(node),
        );
        let (code, message) = match symbol.kind {
            "generic" => {
                if symbol.scope == id {
                    check_generic(project, symbol_id, diagnostics);
                }
                (
                    "OOP007",
                    "public module generic requires a documented exception",
                )
            }
            "procedure" if !symbol.attrs.contains("pointer") => {
                if bound.contains(&symbol_id) {
                    (
                        "OOP006",
                        "type-bound implementation must not be exported directly",
                    )
                } else {
                    (
                        "OOP005",
                        "public non-type-bound procedure requires a documented exception",
                    )
                }
            }
            _ => continue,
        };
        diagnostics.push(
            Diagnostic::new(code, message, location)
                .with_subject(entity_subject(&scope.qualified, name)),
        );
    }
}

fn entity_subject(prefix: &str, entity: &str) -> String {
    if prefix.is_empty() {
        entity.into()
    } else if entity.is_empty() {
        prefix.into()
    } else {
        format!("{prefix}::{entity}")
    }
}

fn context(
    project: &Project<'_>,
    types: &BTreeMap<(usize, usize), &String>,
    unit: usize,
    node: Node<'_>,
) -> String {
    let scope = project
        .scope_of(unit, node)
        .map(|id| &project.scopes[project.scopes[id].declaration.unwrap_or(id)]);
    let mut ancestor = Some(node);
    while let Some(item) = ancestor {
        if let Some(name) = types.get(&(unit, item.id())) {
            return (*name).clone();
        }
        if matches!(
            item.kind(),
            "module" | "submodule" | "program" | "function" | "subroutine" | "module_procedure"
        ) {
            break;
        }
        ancestor = item.parent();
    }
    scope
        .map(|scope| scope.qualified.clone())
        .unwrap_or_default()
}

fn check_polymorphism(
    project: &Project<'_>,
    types: &BTreeMap<(usize, usize), &String>,
    unit: usize,
    node: Node<'_>,
    diagnostics: &mut Vec<Diagnostic>,
) {
    let Some(declaration) = node.parent() else {
        return;
    };
    let prefix = context(project, types, unit, node);
    let source = &project.units[unit];
    let mut emit = |node, name: &str| {
        let subject = entity_subject(&prefix, name);
        diagnostics.push(
            Diagnostic::new(
                "OOP008",
                "CLASS(*) requires a documented exception",
                source.location(node),
            )
            .with_subject(subject),
        );
    };
    if declaration.kind() == "function_statement"
        && let Some(scope) = project
            .scope_of(unit, declaration)
            .map(|id| &project.scopes[id])
    {
        emit(node, scope.result.as_deref().unwrap_or(&scope.name));
        return;
    }
    let names: Vec<_> = fields(declaration, "declarator")
        .into_iter()
        .filter_map(entity)
        .collect();
    if names.is_empty() {
        emit(node, "");
    } else {
        for name in names {
            emit(name, &source.name(name));
        }
    }
}
