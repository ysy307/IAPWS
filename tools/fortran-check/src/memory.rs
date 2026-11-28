use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use tree_sitter::Node;

use crate::{
    config::Config,
    diagnostic::Diagnostic,
    flow::{self, Edge, Step, unwrap},
    library_calls::{
        ArgumentSelector, Call, Storage, arguments, assigned, base, declarator, escapes,
        executable, int, read_call, resource_access, status_tests, storage,
    },
    semantic::{Project, Symbol},
    source::{Unit, children, fields, walk},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Association {
    Undefined,
    Null,
    Unknown,
    Dangling,
    Target(usize),
    Allocated(usize),
    Borrowed(usize),
}

impl Association {
    fn target(self) -> bool {
        matches!(
            self,
            Self::Target(_) | Self::Allocated(_) | Self::Borrowed(_)
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Alternative {
    value: Association,
    guard: Option<(Storage, bool)>,
}

type Facts = BTreeMap<usize, BTreeSet<Alternative>>;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Action {
    None,
    Allocate,
    Deallocate,
    Nullify,
    Associate,
    Borrow,
    Restore,
    Output,
}

struct Operation<'a> {
    action: Action,
    targets: Vec<Node<'a>>,
    status: Option<Storage>,
}

fn symbol<'a>(project: &'_ Project<'a>, scope: usize, node: Node<'a>) -> Option<usize> {
    let unit = project.scopes[scope].unit;
    project
        .expression(unit, node)
        .or_else(|| project.expression(unit, base(node)))
        .or_else(|| storage(project, scope, node).map(|place| place.root))
}

fn attributes(symbol: &Symbol<'_>, names: &[&str]) -> bool {
    names.iter().any(|name| symbol.attrs.contains(*name))
}

fn initial(project: &Project<'_>, scope: usize, pointer: usize) -> BTreeSet<Alternative> {
    let pointer = &project.symbols[pointer];
    let declaration = declarator(pointer);
    let static_local = project.scopes[scope].kind != "program"
        && (pointer.attrs.contains("save") || declaration.kind() == "pointer_init_declarator");
    let values: Vec<_> = if declaration.kind() != "pointer_init_declarator" {
        vec![if static_local {
            Association::Unknown
        } else {
            Association::Undefined
        }]
    } else if children(declaration)
        .last()
        .is_some_and(|node| node.kind() == "null_literal")
    {
        if static_local {
            vec![Association::Null, Association::Unknown]
        } else {
            vec![Association::Null]
        }
    } else {
        vec![Association::Unknown]
    };
    values
        .into_iter()
        .map(|value| Alternative { value, guard: None })
        .collect()
}

fn join(left: &Facts, right: &Facts) -> Facts {
    left.iter()
        .map(|(pointer, values)| (*pointer, values.union(&right[pointer]).cloned().collect()))
        .collect()
}

fn set(values: &BTreeSet<Alternative>, value: Association) -> BTreeSet<Alternative> {
    values
        .iter()
        .map(|old| Alternative {
            value,
            guard: old.guard.clone(),
        })
        .collect()
}

fn guards(facts: &mut Facts, status: &Storage, success: bool) -> bool {
    for alternatives in facts.values_mut() {
        *alternatives = alternatives
            .iter()
            .filter_map(|alternative| {
                if alternative
                    .guard
                    .as_ref()
                    .is_some_and(|(key, sense)| key == status && *sense != success)
                {
                    return None;
                }
                Some(Alternative {
                    value: alternative.value,
                    guard: Some((status.clone(), success)),
                })
            })
            .collect();
    }
    facts.values().all(|values| !values.is_empty())
}

fn is_associated(project: &Project<'_>, scope: usize, node: Node<'_>) -> bool {
    node.kind() == "call_expression"
        && project
            .expression(project.scopes[scope].unit, node)
            .is_some_and(|id| {
                project.symbols[id].kind == "intrinsic" && project.symbols[id].name == "associated"
            })
}

pub(crate) fn formal_arguments<'a>(
    project: &Project<'a>,
    scope: usize,
    node: Node<'a>,
) -> Vec<(Node<'a>, usize)> {
    let item = &project.scopes[scope];
    let unit = &project.units[item.unit];
    let Some(mut id) = project.expression(item.unit, node) else {
        return Vec::new();
    };
    let binding = (project.symbols[id].kind == "binding").then_some(id);
    if let Some(binding) = binding {
        let Some(target) = project.binding_target(binding) else {
            return Vec::new();
        };
        id = target;
    }
    let Some(body) = project.symbols[id]
        .procedure
        .or_else(|| project.implementation(id))
    else {
        return Vec::new();
    };
    let body = &project.scopes[body];
    let signature = body.declaration.map_or(body, |decl| &project.scopes[decl]);
    let mut names = signature.parameters.clone();
    if let Some(binding) = binding.filter(|id| !project.symbols[*id].attrs.contains("nopass")) {
        let passed = project.symbols[binding]
            .attrs
            .iter()
            .find_map(|attr| attr.strip_prefix("pass(")?.strip_suffix(')'));
        if let Some(passed) = passed {
            names.retain(|name| name != passed);
        } else if !names.is_empty() {
            names.remove(0);
        }
    }
    arguments(unit, node)
        .into_iter()
        .enumerate()
        .filter_map(|(index, argument)| {
            let name = argument.keyword.as_ref().or_else(|| names.get(index))?;
            let formal = *signature.symbols.get(name)?;
            Some((argument.node, formal))
        })
        .collect()
}

fn pointer_output<'a>(
    project: &Project<'a>,
    config: &Config,
    scope: usize,
    node: Node<'a>,
) -> Option<Node<'a>> {
    let info = &project.scopes[scope];
    let id = project.expression(info.unit, node)?;
    let item = &project.symbols[id];
    let qualified = item.external.as_ref().map_or_else(
        || item.qualified.clone(),
        |module| format!("{module}::{}", item.name),
    );
    let selector = config.safety.pointer_outputs.get(&qualified).or_else(|| {
        if item.name == "c_f_pointer" && item.external.as_deref() != Some("iso_c_binding") {
            return None;
        }
        config.safety.pointer_outputs.get(&item.name)
    })?;
    argument_at(&project.units[info.unit], node, selector)
}

fn argument_at<'a>(unit: &Unit, node: Node<'a>, selector: &ArgumentSelector) -> Option<Node<'a>> {
    Call {
        family: crate::library_calls::Family::Mpi,
        name: String::new(),
        arguments: arguments(unit, node),
    }
    .argument(selector)
}

fn operation<'a>(
    project: &Project<'a>,
    config: &Config,
    scope: usize,
    node: Node<'a>,
    pointers: &BTreeSet<usize>,
) -> Operation<'a> {
    let unit = &project.units[project.scopes[scope].unit];
    let targets = |node| {
        children(node)
            .into_iter()
            .filter(|item| item.kind() != "keyword_argument")
            .collect()
    };
    let status = || {
        children(node).into_iter().find_map(|item| {
            let key = item.child_by_field_name("name")?;
            if item.kind() != "keyword_argument" || unit.name(key) != "stat" {
                return None;
            }
            storage(project, scope, item.child_by_field_name("value")?)
        })
    };
    match node.kind() {
        "allocate_statement" => Operation {
            action: Action::Allocate,
            targets: fields(node, "allocation")
                .into_iter()
                .map(|item| {
                    if item.kind() == "sized_allocation" {
                        item.named_child(0).unwrap_or(item)
                    } else {
                        item
                    }
                })
                .collect(),
            status: status(),
        },
        "deallocate_statement" => Operation {
            action: Action::Deallocate,
            targets: targets(node),
            status: status(),
        },
        "nullify_statement" => Operation {
            action: Action::Nullify,
            targets: targets(node),
            status: None,
        },
        "pointer_association_statement" => Operation {
            action: Action::Associate,
            targets: targets(node),
            status: None,
        },
        "subroutine_call" => {
            if let Some(target) = pointer_output(project, config, scope, node) {
                return Operation {
                    action: Action::Output,
                    targets: vec![target],
                    status: None,
                };
            }
            if let Some(call) = read_call(project, config, scope, node) {
                if let Some(resource) =
                    resource_access(&config.safety, &call).filter(|item| item.borrowed)
                {
                    return Operation {
                        action: if resource.acquire {
                            Action::Borrow
                        } else {
                            Action::Restore
                        },
                        targets: vec![resource.handle],
                        status: resource
                            .status
                            .and_then(|status| storage(project, scope, status)),
                    };
                }
                return Operation {
                    action: Action::None,
                    targets: Vec::new(),
                    status: None,
                };
            }
            let mut targets: Vec<_> = formal_arguments(project, scope, node)
                .into_iter()
                .filter_map(|(actual, id)| {
                    let formal = &project.symbols[id];
                    (formal.attrs.contains("pointer")
                        && attributes(formal, &["intent(out)", "intent(inout)"]))
                    .then_some(actual)
                })
                .collect();
            let opaque = project
                .expression(project.scopes[scope].unit, node)
                .is_none_or(|id| project.symbols[id].kind == "external");
            if opaque {
                targets.extend(arguments(unit, node).into_iter().filter_map(|argument| {
                    let id = symbol(project, scope, argument.node)?;
                    (argument.node.kind() == "identifier" && pointers.contains(&id))
                        .then_some(argument.node)
                }));
            }
            Operation {
                action: Action::Output,
                targets,
                status: None,
            }
        }
        _ => Operation {
            action: Action::None,
            targets: Vec::new(),
            status: None,
        },
    }
}

fn release(facts: &mut Facts, holder: usize, restoring: bool) {
    let Some(values) = facts.get(&holder) else {
        return;
    };
    let candidates: BTreeSet<_> = values.iter().map(|item| item.value).collect();
    for (pointer, values) in facts.iter_mut() {
        if *pointer == holder {
            *values = set(
                values,
                if restoring {
                    Association::Dangling
                } else {
                    Association::Null
                },
            );
            continue;
        }
        let mut next = BTreeSet::new();
        for alternative in values.iter() {
            let invalidated = alternative.value.target() && candidates.contains(&alternative.value);
            if !invalidated || candidates.len() != 1 {
                next.insert(alternative.clone());
            }
            if invalidated {
                next.insert(Alternative {
                    value: Association::Dangling,
                    guard: alternative.guard.clone(),
                });
            }
        }
        *values = next;
    }
}

struct Procedure<'a, 'p> {
    project: &'p Project<'a>,
    config: &'p Config,
    scope: usize,
    unit: &'p Unit,
    pointers: BTreeSet<usize>,
    operations: HashMap<Node<'a>, Operation<'a>>,
    statuses: BTreeSet<Storage>,
}

impl<'a> Procedure<'a, '_> {
    fn location_issue(
        &self,
        code: &str,
        message: String,
        node: Node<'_>,
        pointer: Option<usize>,
        certain: bool,
    ) -> Diagnostic {
        let subject = pointer.map_or(&self.project.scopes[self.scope].qualified, |id| {
            &self.project.symbols[id].qualified
        });
        let diagnostic =
            Diagnostic::new(code, message, self.unit.location(node)).with_subject(subject);
        if certain {
            diagnostic
        } else {
            diagnostic.warning()
        }
    }

    fn symbol(&self, node: Node<'a>) -> Option<usize> {
        symbol(self.project, self.scope, node)
    }

    fn place(&self, node: Node<'_>) -> Option<Storage> {
        storage(self.project, self.scope, node)
    }

    fn transfer(&self, step: &Step<'a>, facts: &Facts) -> Facts {
        let Some(node) = step.node else {
            return facts.clone();
        };
        let operation = &self.operations[&node];
        let mut facts = facts.clone();
        let mut writes = BTreeSet::new();
        if let Some((left, _)) = assigned(node) {
            writes.extend(self.place(left));
        }
        writes.extend(operation.status.clone());
        if node.kind() == "subroutine_call" {
            let resolved = read_call(self.project, self.config, self.scope, node);
            if let Some(call) = resolved {
                if let Some(status) = self
                    .config
                    .safety
                    .expected_status(&call)
                    .and_then(|selector| call.argument(selector))
                {
                    writes.extend(self.place(status));
                }
            } else {
                writes.extend(
                    arguments(self.unit, node)
                        .into_iter()
                        .filter_map(|arg| self.place(arg.node)),
                );
            }
        }
        if !writes.is_disjoint(&self.statuses) {
            for alternatives in facts.values_mut() {
                *alternatives = alternatives
                    .iter()
                    .map(|old| {
                        let guard = old
                            .guard
                            .clone()
                            .filter(|(status, _)| !writes.contains(status));
                        Alternative {
                            value: old.value,
                            guard,
                        }
                    })
                    .collect();
            }
        }
        let mut changed = facts.clone();
        if operation.action == Action::Associate {
            if let [destination, source] = operation.targets.as_slice()
                && let Some(id) = self
                    .symbol(*destination)
                    .filter(|id| self.pointers.contains(id))
            {
                let origin = self.symbol(*source);
                let alternatives = if source.kind() == "null_literal" {
                    set(&facts[&id], Association::Null)
                } else if let Some(values) = origin.and_then(|id| facts.get(&id)) {
                    values.clone()
                } else {
                    set(
                        &facts[&id],
                        origin
                            .filter(|id| {
                                self.project.symbols[*id].kind == "variable"
                                    && !self.project.symbols[*id].attrs.contains("pointer")
                            })
                            .map_or(Association::Unknown, Association::Target),
                    )
                };
                changed.insert(id, alternatives);
            }
            return changed;
        }
        for target in &operation.targets {
            let Some(id) = self.symbol(*target).filter(|id| self.pointers.contains(id)) else {
                continue;
            };
            if matches!(operation.action, Action::Restore | Action::Deallocate) {
                release(&mut changed, id, operation.action == Action::Restore);
                continue;
            }
            let value = match operation.action {
                Action::Allocate => Association::Allocated(target.start_byte()),
                Action::Borrow => Association::Borrowed(target.start_byte()),
                Action::Nullify => Association::Null,
                _ => Association::Unknown,
            };
            changed.insert(id, set(&facts[&id], value));
        }
        if let Some(status) = &operation.status {
            guards(&mut changed, status, true);
            guards(&mut facts, status, false);
            return join(&changed, &facts);
        }
        changed
    }

    fn refine(&self, edge: &Edge<'a>, facts: &Facts) -> Option<Facts> {
        let mut facts = facts.clone();
        for (status, success) in status_tests(self.project, self.config, self.scope, edge) {
            if self.statuses.contains(&status) && !guards(&mut facts, &status, success) {
                return None;
            }
        }
        if let Some(condition) = edge.condition {
            let mut call = unwrap(condition);
            let mut truth = edge.truth;
            if self.unit.name(call).starts_with(".not.") {
                call = unwrap(call.named_child(0).unwrap_or(call));
                truth = !truth;
            }
            if is_associated(self.project, self.scope, call) {
                let args = arguments(self.unit, call);
                if let [actual] = args.as_slice()
                    && let Some(values) = self.symbol(actual.node).and_then(|id| facts.get_mut(&id))
                {
                    values.retain(|value| {
                        if truth {
                            value.value != Association::Null
                        } else {
                            !value.value.target()
                        }
                    });
                    if values.is_empty() {
                        return None;
                    }
                }
            }
        }
        Some(facts)
    }

    fn pointer_uses(&self, node: Node<'a>, facts: &Facts, diagnostics: &mut Vec<Diagnostic>) {
        if matches!(
            node.kind(),
            "variable_declaration" | "variable_modification" | "save_statement"
        ) {
            return;
        }
        let operation = &self.operations[&node];
        let mut descriptor = HashSet::new();
        let mut definitions = HashSet::new();
        let mut opaque = HashSet::new();
        if operation.action == Action::Associate {
            if let [destination, source] = operation.targets.as_slice() {
                definitions.insert(base(*destination));
                if matches!(
                    source.kind(),
                    "identifier" | "derived_type_member_expression"
                ) {
                    descriptor.insert(base(*source));
                }
            }
        } else if matches!(
            operation.action,
            Action::Allocate | Action::Nullify | Action::Borrow | Action::Output | Action::Restore
        ) {
            definitions.extend(operation.targets.iter().map(|target| base(*target)));
        }
        if node.kind() == "subroutine_call" {
            let formals: HashMap<_, _> = formal_arguments(self.project, self.scope, node)
                .into_iter()
                .collect();
            for item in arguments(self.unit, node) {
                if item.node.kind() == "identifier"
                    && formals
                        .get(&item.node)
                        .is_none_or(|id| self.project.symbols[*id].attrs.contains("pointer"))
                {
                    opaque.insert(item.node);
                }
            }
        }
        for item in walk(node) {
            if is_associated(self.project, self.scope, item) {
                descriptor.extend(
                    arguments(self.unit, item)
                        .into_iter()
                        .filter_map(|arg| (arg.node.kind() == "identifier").then_some(arg.node)),
                );
            }
            if item.kind() != "identifier" || definitions.contains(&item) {
                continue;
            }
            let Some(id) = self.symbol(item) else {
                continue;
            };
            let Some(values) = facts.get(&id) else {
                continue;
            };
            let pointer = &self.project.symbols[id];
            let undefined = values
                .iter()
                .filter(|value| value.value == Association::Undefined)
                .count();
            if undefined != 0 && !opaque.contains(&item) {
                diagnostics.push(self.location_issue(
                    "MEM002",
                    format!(
                        "Pointer '{}' may have undefined association status",
                        pointer.name
                    ),
                    item,
                    Some(id),
                    undefined == values.len(),
                ));
            }
            let invalid = values
                .iter()
                .filter(|value| {
                    value.value == Association::Dangling
                        || (value.value == Association::Null
                            && !descriptor.contains(&item)
                            && !opaque.contains(&item))
                })
                .count();
            if invalid != 0 {
                diagnostics.push(self.location_issue(
                    "MEM003",
                    format!(
                        "Pointer '{}' may reference released or absent storage",
                        pointer.name
                    ),
                    item,
                    Some(id),
                    invalid == values.len(),
                ));
            }
            if operation.action == Action::Deallocate && operation.targets.contains(&item) {
                let static_targets = values
                    .iter()
                    .filter(|value| {
                        let Association::Target(target) = value.value else {
                            return false;
                        };
                        let target = &self.project.symbols[target];
                        target.scope == self.scope
                            && !target.attrs.contains("allocatable")
                            && !self.project.scopes[self.scope]
                                .parameters
                                .contains(&target.name)
                    })
                    .count();
                if static_targets != 0 {
                    diagnostics.push(self.location_issue(
                        "MEM003",
                        format!(
                            "Pointer '{}' may deallocate nonallocated TARGET storage",
                            pointer.name
                        ),
                        item,
                        Some(id),
                        static_targets == values.len(),
                    ));
                }
            }
        }
    }

    fn other_checks(&self, node: Node<'a>, diagnostics: &mut Vec<Diagnostic>) {
        if node.kind() == "assignment_statement"
            && let Some((left, right)) = assigned(node)
        {
            let symbol = self.symbol(left).map(|id| &self.project.symbols[id]);
            if let Some(symbol) = symbol.filter(|symbol| {
                left.kind() != "call_expression" && symbol.attrs.contains("allocatable")
            }) {
                let declaration = declarator(symbol).parent();
                let dtype = declaration.and_then(|node| node.child_by_field_name("type"));
                let deferred = dtype
                    .map(|node| {
                        self.project.units[self.project.scopes[symbol.scope].unit].name(node)
                    })
                    .is_some_and(|name| {
                        name.starts_with("character(len=:") || name.starts_with("character(:")
                    });
                let scalar = matches!(
                    right.kind(),
                    "number_literal" | "boolean_literal" | "string_literal"
                ) || (right.kind() == "identifier"
                    && self
                        .symbol(right)
                        .is_some_and(|id| !ranked(&self.project.symbols[id])));
                if (deferred || (ranked(symbol) && !scalar))
                    && self.unit.name(left) != self.unit.name(right)
                {
                    diagnostics.push(self.location_issue(
                        "MEM004",
                        format!(
                            "Whole allocatable assignment to '{}' may reallocate storage",
                            symbol.name
                        ),
                        left,
                        None,
                        false,
                    ));
                }
            }
        }
        for call in
            walk(node).filter(|node| matches!(node.kind(), "subroutine_call" | "call_expression"))
        {
            for (actual, formal) in formal_arguments(self.project, self.scope, call) {
                let formal = &self.project.symbols[formal];
                let formal_unit = &self.project.units[self.project.scopes[formal.scope].unit];
                let Some(id) = self.symbol(actual) else {
                    continue;
                };
                let symbol = &self.project.symbols[id];
                if actual.kind() == "call_expression"
                    && matches!(symbol.kind, "variable" | "constant")
                    && packed_formal(formal_unit, formal)
                    && strided(self.project, self.unit, actual, symbol)
                {
                    diagnostics.push(self.location_issue(
                        "MEM005",
                        format!(
                            "Strided actual argument for '{}' may require an array temporary",
                            formal.name
                        ),
                        actual,
                        None,
                        false,
                    ));
                }
            }
        }
    }
}

fn ranked(symbol: &Symbol<'_>) -> bool {
    symbol
        .attrs
        .iter()
        .any(|attr| attr.starts_with("dimension("))
        || walk(declarator(symbol)).any(|node| node.kind() == "sized_declarator")
}

fn packed_formal(unit: &Unit, symbol: &Symbol<'_>) -> bool {
    if attributes(symbol, &["pointer", "allocatable"]) {
        return false;
    }
    let size = walk(declarator(symbol)).find(|node| node.kind() == "size");
    if let Some(size) = size {
        let bounds = children(size);
        if bounds
            .iter()
            .any(|item| item.kind() == "assumed_rank" || unit.name(*item) == "..")
        {
            return false;
        }
        return symbol.attrs.contains("contiguous")
            || bounds.iter().all(|item| !unit.name(*item).ends_with(':'));
    }
    ranked(symbol) && symbol.attrs.contains("contiguous")
}

fn strided(project: &Project<'_>, unit: &Unit, actual: Node<'_>, symbol: &Symbol<'_>) -> bool {
    let declared_unit = &project.units[project.scopes[symbol.scope].unit];
    let bounds = walk(declarator(symbol))
        .find(|node| node.kind() == "size")
        .map(children)
        .unwrap_or_default();
    let mut gaps = false;
    let mut strided = false;
    for (index, subscript) in arguments(unit, actual).into_iter().enumerate() {
        let text = unit.name(subscript.node);
        let declared = bounds
            .get(index)
            .map_or_else(String::new, |bound| declared_unit.name(*bound));
        let (lo, hi) = if let Some((low, high)) = declared.split_once(':') {
            (if low.is_empty() { Some(1) } else { int(low) }, int(high))
        } else {
            (Some(1), int(&declared))
        };
        let extent = lo
            .zip(hi)
            .and_then(|(lo, hi)| hi.checked_sub(lo)?.checked_add(1));
        if subscript.node.kind() != "extent_specifier" {
            gaps |= extent.is_some_and(|value| value > 1) && int(&text).is_some();
            continue;
        }
        let parts: Vec<_> = text.split(':').collect();
        let [first, last, ..] = parts.as_slice() else {
            continue;
        };
        let Some(step) = parts.get(2).map_or(Some(1), |step| int(step)) else {
            continue;
        };
        if step == 0 {
            continue;
        }
        let first = if first.is_empty() { lo } else { int(first) };
        let last = if last.is_empty() { hi } else { int(last) };
        let count = first.zip(last).and_then(|(first, last)| {
            if (step > 0 && last < first) || (step < 0 && last > first) {
                return Some(0);
            }
            last.checked_sub(first)?.checked_div(step)?.checked_add(1)
        });
        if count == Some(0) {
            return false;
        }
        strided |= (step != 1 || gaps) && count != Some(1);
        gaps |= extent
            .zip(count)
            .is_some_and(|(extent, count)| 0 < count && count < extent);
    }
    strided
}

pub fn check(project: &Project<'_>, config: &Config) -> Vec<Diagnostic> {
    let mut diagnostics = Vec::new();
    for (scope, info) in project
        .scopes
        .iter()
        .enumerate()
        .filter(|(_, info)| executable(info))
    {
        let unit = &project.units[info.unit];
        let pointers: BTreeSet<_> = info
            .symbols
            .values()
            .copied()
            .filter(|id| {
                let symbol = &project.symbols[*id];
                symbol.scope == scope
                    && symbol.kind == "variable"
                    && symbol.attrs.contains("pointer")
                    && !info.parameters.contains(&symbol.name)
                    && info.result.as_ref() != Some(&symbol.name)
            })
            .collect();
        let flow = flow::build(info, unit, |node| {
            config
                .safety
                .terminates(read_call(project, config, scope, node).as_ref())
        });
        let operations: HashMap<_, _> = flow
            .steps
            .iter()
            .filter_map(|step| {
                let node = step.node?;
                Some((node, operation(project, config, scope, node, &pointers)))
            })
            .collect();
        let statuses = operations
            .values()
            .filter_map(|operation| operation.status.clone())
            .collect();
        let procedure = Procedure {
            project,
            config,
            scope,
            unit,
            pointers,
            operations,
            statuses,
        };
        let initial = procedure
            .pointers
            .iter()
            .map(|pointer| (*pointer, initial(project, scope, *pointer)))
            .collect();
        let facts = flow::solve(
            &flow,
            initial,
            |step, facts| procedure.transfer(step, facts),
            join,
            |edge, facts| procedure.refine(edge, facts),
        );
        let mut allocated = BTreeSet::new();
        let mut descriptors = BTreeSet::new();
        for (step, facts) in flow.steps.iter().zip(facts) {
            let (Some(node), Some(facts)) = (step.node, facts) else {
                continue;
            };
            procedure.pointer_uses(node, &facts, &mut diagnostics);
            procedure.other_checks(node, &mut diagnostics);
            let operation = &procedure.operations[&node];
            if operation.action == Action::Allocate {
                allocated.extend(
                    operation
                        .targets
                        .iter()
                        .filter_map(|target| procedure.symbol(*target)),
                );
            } else if matches!(
                operation.action,
                Action::Associate
                    | Action::Nullify
                    | Action::Output
                    | Action::Restore
                    | Action::Borrow
            ) {
                descriptors.extend(
                    operation
                        .targets
                        .iter()
                        .filter_map(|target| procedure.symbol(base(*target))),
                );
            }
            for call in walk(node)
                .filter(|node| matches!(node.kind(), "subroutine_call" | "call_expression"))
            {
                if project
                    .expression(info.unit, call)
                    .is_some_and(|id| project.symbols[id].kind != "variable")
                {
                    descriptors.extend(
                        arguments(unit, call)
                            .into_iter()
                            .filter_map(|arg| procedure.symbol(base(arg.node))),
                    );
                }
            }
        }
        for id in procedure
            .pointers
            .intersection(&allocated)
            .filter(|id| !descriptors.contains(id))
        {
            let pointer = &project.symbols[*id];
            if !escapes(
                project,
                scope,
                &Storage {
                    root: *id,
                    parts: Vec::new(),
                },
            ) {
                diagnostics.push(procedure.location_issue(
                    "MEM001",
                    format!(
                        "Locally owned pointer '{}' can use ALLOCATABLE",
                        pointer.name
                    ),
                    pointer.node,
                    Some(*id),
                    false,
                ));
            }
        }
    }
    diagnostics
}
