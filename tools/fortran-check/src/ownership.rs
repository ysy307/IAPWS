use std::collections::{BTreeMap, BTreeSet, HashMap};

use tree_sitter::Node;

use crate::{
    config::Config,
    diagnostic::{Diagnostic, Location},
    flow::{self, Edge, Flow, Step},
    library_calls::{
        Part, ResourceAccess, Storage, assigned, base, body_nodes, declarator, escapes, executable,
        identity, read_call, resource_access, status_tests, storage,
    },
    semantic::Project,
    source::{Unit, children, fields},
};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Live {
    handle: Storage,
    release: String,
    context: Vec<Option<String>>,
    borrowed: bool,
    site: Location,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Pending {
    Acquire(Live),
    Release(Live, bool),
}

impl Pending {
    fn resource(&self) -> &Live {
        match self {
            Self::Acquire(resource) | Self::Release(resource, _) => resource,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct State {
    live: BTreeSet<Live>,
    definite: BTreeSet<Live>,
    pending: BTreeSet<(Storage, Pending)>,
}

impl State {
    fn join(left: &Self, right: &Self) -> Self {
        Self {
            live: left.live.union(&right.live).cloned().collect(),
            definite: left
                .definite
                .intersection(&right.definite)
                .cloned()
                .collect(),
            pending: left.pending.intersection(&right.pending).cloned().collect(),
        }
    }

    fn forget(&mut self, lost: &BTreeSet<Live>) {
        self.live.retain(|item| !lost.contains(item));
        self.definite.retain(|item| !lost.contains(item));
        self.pending
            .retain(|(_, pending)| !lost.contains(pending.resource()));
    }

    fn overwritten(&mut self, storage: &Storage) {
        self.pending.retain(|(key, _)| key != storage);
    }
}

struct Procedure<'a, 'p> {
    project: &'p Project<'a>,
    scope: usize,
    config: &'p Config,
    unit: &'p Unit,
    flow: Flow<'a>,
    calls: HashMap<Node<'a>, (Option<Storage>, Option<ResourceAccess<'a>>)>,
}

impl<'a, 'p> Procedure<'a, 'p> {
    fn new(project: &'p Project<'a>, config: &'p Config, scope: usize) -> Self {
        let info = &project.scopes[scope];
        let unit = &project.units[info.unit];
        let calls: HashMap<_, _> = body_nodes(info)
            .into_iter()
            .filter_map(|node| read_call(project, config, scope, node).map(|call| (node, call)))
            .collect();
        let flow = flow::build(info, unit, |node| {
            config.safety.terminates(calls.get(&node))
        });
        let accesses = calls
            .into_iter()
            .map(|(node, call)| {
                let resource = resource_access(&config.safety, &call);
                let status = resource.as_ref().and_then(|item| item.status).or_else(|| {
                    config
                        .safety
                        .expected_status(&call)
                        .and_then(|selector| call.argument(selector))
                });
                let status = status.and_then(|status| storage(project, scope, status));
                (node, (status, resource))
            })
            .collect();
        Self {
            project,
            scope,
            config,
            unit,
            flow,
            calls: accesses,
        }
    }

    fn place(&self, node: Node<'_>) -> Option<Storage> {
        storage(self.project, self.scope, node).filter(|place| {
            !place
                .parts
                .iter()
                .any(|part| matches!(part, Part::Indices(items) if items.contains(&None)))
        })
    }

    fn identity(&self, operation: &ResourceAccess<'a>) -> Vec<Option<String>> {
        operation
            .context
            .iter()
            .map(|node| node.map(|node| identity(self.project, self.scope, node)))
            .collect()
    }

    fn loc(&self, code: &str, message: String, location: Location, certain: bool) -> Diagnostic {
        let diagnostic = Diagnostic::new(code, message, location)
            .with_subject(&self.project.scopes[self.scope].qualified);
        if certain {
            diagnostic
        } else {
            diagnostic.warning()
        }
    }

    fn diagnostic(&self, code: &str, message: String, node: Node<'_>, certain: bool) -> Diagnostic {
        self.loc(code, message, self.unit.location(node), certain)
    }

    fn type_at(&self, node: Node<'_>) -> Option<usize> {
        self.project
            .expression(self.project.scopes[self.scope].unit, node)
            .or_else(|| self.place(node).map(|place| place.root))
            .and_then(|id| self.project.type_of(id))
    }

    fn defined_assignment(&self, left: Node<'_>, right: Node<'_>) -> bool {
        self.project
            .lookup_kind(self.scope, "assignment(=)", &["generic"])
            .is_some()
            || [left, right]
                .into_iter()
                .filter_map(|node| self.type_at(node))
                .any(|owner| self.project.member(owner, "assignment(=)").is_some())
    }

    fn overwrite(&self, state: &State, storage: &Storage, left: Node<'_>) -> BTreeSet<Live> {
        let pointer_left = self
            .project
            .expression(self.project.scopes[self.scope].unit, left)
            .is_some_and(|id| self.project.symbols[id].attrs.contains("pointer"));
        state
            .live
            .iter()
            .filter(|resource| &resource.handle == storage && (!resource.borrowed || !pointer_left))
            .cloned()
            .collect()
    }

    fn transfer(&self, step: &Step<'a>, before: &State) -> State {
        let mut after = before.clone();
        let Some(node) = step.node else {
            return after;
        };
        if let Some((status, operation)) = self.calls.get(&node) {
            if let Some(status) = status {
                after.overwritten(status);
            }
            let Some(operation) = operation else {
                return after;
            };
            let Some(handle) = self.place(operation.handle) else {
                return after;
            };
            let context = self.identity(operation);
            let affected: BTreeSet<_> = before
                .live
                .iter()
                .filter(|resource| {
                    resource.handle == handle
                        && (operation.acquire || matching(resource, operation, &context))
                })
                .cloned()
                .collect();
            after.forget(&affected);
            if operation.acquire {
                let resource = Live {
                    handle,
                    release: operation.release.clone(),
                    context,
                    borrowed: operation.borrowed,
                    site: self.unit.location(operation.handle),
                };
                after.live.insert(resource.clone());
                if let Some(status) = status {
                    after
                        .pending
                        .insert((status.clone(), Pending::Acquire(resource)));
                } else {
                    after.definite.insert(resource);
                }
            }
            if let Some(status) = status {
                after.pending.extend(affected.into_iter().map(|resource| {
                    let was_definite = before.definite.contains(&resource);
                    (status.clone(), Pending::Release(resource, was_definite))
                }));
            }
            return after;
        }
        let Some((left, right)) = assigned(node).filter(|_| node.kind() == "assignment_statement")
        else {
            return after;
        };
        let Some(destination) = self.place(left) else {
            return after;
        };
        let source = self.place(right);
        after.overwritten(&destination);
        if source.as_ref() == Some(&destination) || self.defined_assignment(left, right) {
            return after;
        }
        after.forget(&self.overwrite(before, &destination, left));
        if let Some(source) = source.filter(|_| escapes(self.project, self.scope, &destination)) {
            let escaped = after
                .live
                .iter()
                .filter(|resource| resource.handle == source && !resource.borrowed)
                .cloned()
                .collect();
            after.forget(&escaped);
        }
        after
    }

    fn refine(&self, edge: &Edge<'a>, state: &State) -> Option<State> {
        let mut state = state.clone();
        for (status, succeeded) in status_tests(self.project, self.config, self.scope, edge) {
            let consequences: Vec<_> = state
                .pending
                .iter()
                .filter(|(place, _)| *place == status)
                .map(|(_, item)| item.clone())
                .collect();
            state.overwritten(&status);
            for pending in consequences {
                match (pending, succeeded) {
                    (Pending::Acquire(item), true) => {
                        state.definite.insert(item);
                    }
                    (Pending::Acquire(item), false) => {
                        state.live.remove(&item);
                        state.definite.remove(&item);
                    }
                    (Pending::Release(item, definite), false) => {
                        state.live.insert(item.clone());
                        if definite {
                            state.definite.insert(item);
                        }
                    }
                    (Pending::Release(_, _), true) => {}
                }
            }
        }
        Some(state)
    }

    fn diagnostics(&self, owners: &BTreeSet<usize>, diagnostics: &mut Vec<Diagnostic>) {
        let facts = flow::solve(
            &self.flow,
            State::default(),
            |step, before| self.transfer(step, before),
            State::join,
            |edge, state| self.refine(edge, state),
        );
        if let Some(facts) = &facts[self.flow.exit] {
            for resource in &facts.live {
                if !escapes(self.project, self.scope, &resource.handle) {
                    let handle = &self.project.symbols[resource.handle.root].name;
                    diagnostics.push(self.loc(
                        "OWN001",
                        format!(
                            "Resource in '{handle}' may reach a normal exit without {}",
                            resource.release.to_ascii_uppercase()
                        ),
                        resource.site.clone(),
                        false,
                    ));
                }
            }
        }
        for (step, before) in self.flow.steps.iter().zip(facts) {
            let (Some(node), Some(before)) = (step.node, before) else {
                continue;
            };
            if let Some((_, Some(operation))) = self.calls.get(&node) {
                self.check_operation(operation, &before, diagnostics);
            } else if node.kind() == "assignment_statement" {
                self.check_assignment(node, &before, owners, diagnostics);
            }
        }
    }

    fn check_operation(
        &self,
        operation: &ResourceAccess<'a>,
        before: &State,
        diagnostics: &mut Vec<Diagnostic>,
    ) {
        let Some(handle) = self.place(operation.handle) else {
            return;
        };
        let active: Vec<_> = before
            .live
            .iter()
            .filter(|resource| resource.handle == handle)
            .collect();
        if active.is_empty() {
            return;
        }
        let name = self.unit.name(operation.handle);
        if operation.acquire {
            diagnostics.push(self.diagnostic(
                "OWN005",
                format!("Acquisition can overwrite the active resource in '{name}'"),
                operation.handle,
                false,
            ));
            return;
        }
        let context = self.identity(operation);
        let wrong: Vec<_> = active
            .into_iter()
            .filter(|resource| !matching(resource, operation, &context))
            .collect();
        if wrong.is_empty() {
            return;
        }
        let certain = wrong.iter().any(|resource| {
            before.definite.contains(*resource) && resource.release != operation.release
        });
        let expected: BTreeSet<_> = wrong
            .iter()
            .map(|resource| resource.release.to_ascii_uppercase())
            .collect();
        let expected = expected.into_iter().collect::<Vec<_>>().join(", ");
        diagnostics.push(self.diagnostic(
            "OWN002",
            format!(
                "Release of '{name}' does not match the acquisition and its context ({expected})"
            ),
            operation.handle,
            certain,
        ));
    }

    fn check_assignment(
        &self,
        node: Node<'a>,
        before: &State,
        owners: &BTreeSet<usize>,
        diagnostics: &mut Vec<Diagnostic>,
    ) {
        let Some((left, right)) = assigned(node) else {
            return;
        };
        let Some(destination) = self.place(left) else {
            return;
        };
        let source = self.place(right);
        if source.as_ref() == Some(&destination) || self.defined_assignment(left, right) {
            return;
        }
        if !self.overwrite(before, &destination, left).is_empty() {
            let name = self.unit.name(left);
            diagnostics.push(self.diagnostic(
                "OWN005",
                format!("Assignment can overwrite the active resource in '{name}'"),
                left,
                false,
            ));
        }
        let Some(source) = source else {
            return;
        };
        let copied = self.type_at(right).is_some_and(|id| owners.contains(&id));
        let active = before.live.iter().any(|resource| resource.handle == source)
            && !escapes(self.project, self.scope, &destination);
        if copied || active {
            diagnostics.push(self.diagnostic(
                "OWN003",
                "Intrinsic assignment can copy an owning value or acquired handle".into(),
                node,
                false,
            ));
        }
    }
}

fn matching(resource: &Live, operation: &ResourceAccess<'_>, context: &[Option<String>]) -> bool {
    resource.release == operation.release
        && resource.context.len() == context.len()
        && resource
            .context
            .iter()
            .zip(context)
            .all(|(left, right)| left.is_none() || right.is_none() || left == right)
}

type FieldsByOwner = BTreeMap<usize, BTreeSet<usize>>;

fn field_events(
    project: &Project<'_>,
    procedures: &[Procedure<'_, '_>],
) -> (FieldsByOwner, FieldsByOwner) {
    let mut owned = FieldsByOwner::new();
    let mut released = FieldsByOwner::new();
    for procedure in procedures {
        let scope = &project.scopes[procedure.scope];
        for node in procedure.flow.steps.iter().filter_map(|step| step.node) {
            let mut fields = Vec::new();
            let mut pointer = false;
            if let Some((_, Some(resource))) = procedure.calls.get(&node).filter(|(_, resource)| {
                resource.as_ref().is_some_and(|resource| !resource.borrowed)
            }) {
                fields.push((resource.handle, resource.acquire));
            } else if node.kind() == "allocate_statement" || node.kind() == "deallocate_statement" {
                pointer = true;
                if node.kind() == "allocate_statement" {
                    fields.extend(fields_for_alloc(node).into_iter().map(|node| (node, true)));
                } else {
                    fields.extend(
                        children(node)
                            .into_iter()
                            .filter(|node| node.kind() != "keyword_argument")
                            .map(|node| (node, false)),
                    );
                }
            }
            for (handle, acquires) in fields {
                let Some(field) = project.expression(scope.unit, handle) else {
                    continue;
                };
                let symbol = &project.symbols[field];
                let Some(owner) = symbol.owner else {
                    continue;
                };
                if pointer && !symbol.attrs.contains("pointer") {
                    continue;
                }
                if acquires {
                    owned.entry(owner).or_default().insert(field);
                } else if let Some(receiver) = project.expression(scope.unit, base(handle))
                    && scope.parameters.contains(&project.symbols[receiver].name)
                {
                    released.entry(procedure.scope).or_default().insert(field);
                }
            }
        }
    }
    (owned, released)
}

fn fields_for_alloc(node: Node<'_>) -> Vec<Node<'_>> {
    fields(node, "allocation")
        .into_iter()
        .map(|node| {
            if node.kind() == "sized_allocation" {
                node.named_child(0).unwrap_or(node)
            } else {
                node
            }
        })
        .collect()
}

fn owners(project: &Project<'_>, fields: &FieldsByOwner) -> BTreeSet<usize> {
    let mut result: BTreeSet<_> = fields.keys().copied().collect();
    loop {
        let added: Vec<_> = project
            .types
            .iter()
            .copied()
            .filter(|owner| {
                !result.contains(owner)
                    && (project
                        .resolve_parent(*owner)
                        .is_some_and(|parent| result.contains(&parent))
                        || project.symbols[*owner].members.values().any(|field| {
                            let symbol = &project.symbols[*field];
                            symbol.kind == "variable"
                                && !symbol.attrs.contains("pointer")
                                && project
                                    .type_of(*field)
                                    .is_some_and(|dtype| result.contains(&dtype))
                        }))
            })
            .collect();
        if added.is_empty() {
            return result;
        }
        result.extend(added);
    }
}

fn cleanups(
    project: &Project<'_>,
    owned: &FieldsByOwner,
    released: &FieldsByOwner,
    diagnostics: &mut Vec<Diagnostic>,
) {
    for owner in &project.types {
        let symbol = &project.symbols[*owner];
        let unit = &project.units[project.scopes[symbol.scope].unit];
        let mut cleaned = BTreeSet::new();
        for binding in symbol.members.values().copied() {
            let item = &project.symbols[binding];
            if !matches!(item.kind, "binding" | "final") {
                continue;
            }
            let Some(target) = project.binding_target(binding) else {
                continue;
            };
            let Some(scope) = project.implementation(target) else {
                continue;
            };
            let body = &project.scopes[scope];
            if item.kind == "final" {
                let signature = body
                    .declaration
                    .map_or(body, |declaration| &project.scopes[declaration]);
                let valid = signature.parameters.len() == 1
                    && matches!(body.kind, "subroutine" | "module_procedure")
                    && signature
                        .symbols
                        .get(&signature.parameters[0])
                        .is_some_and(|id| {
                            let dummy = &project.symbols[*id];
                            let type_decl = declarator(dummy)
                                .parent()
                                .and_then(|node| node.child_by_field_name("type"));
                            type_decl.is_some_and(|node| {
                                project.units[signature.unit]
                                    .name(node)
                                    .starts_with("type(")
                            }) && project.type_of(*id) == Some(*owner)
                                && !["optional", "pointer", "allocatable", "value", "intent(out)"]
                                    .iter()
                                    .any(|attr| dummy.attrs.contains(*attr))
                        });
                if !valid {
                    diagnostics.push(
                        Diagnostic::new(
                            "OWN004",
                            format!("Invalid FINAL binding for '{}'", symbol.name),
                            unit.location(item.node),
                        )
                        .with_subject(&item.qualified),
                    );
                    continue;
                }
            }
            if let Some(fields) = released.get(&scope) {
                cleaned.extend(fields.iter().copied());
            }
        }
        if let Some(fields) = owned.get(owner) {
            for field in fields.difference(&cleaned) {
                let field = &project.symbols[*field];
                diagnostics.push(Diagnostic::new("OWN004", format!("Owned component '{}' has no recognized cleanup in a FINAL or type-bound procedure", field.name),
                    unit.location(field.node)).with_subject(&field.qualified).warning());
            }
        }
    }
}

pub fn check(project: &Project<'_>, config: &Config) -> Vec<Diagnostic> {
    let procedures: Vec<_> = project
        .scopes
        .iter()
        .enumerate()
        .filter(|(_, scope)| executable(scope))
        .map(|(scope, _)| Procedure::new(project, config, scope))
        .collect();
    let (owned, released) = field_events(project, &procedures);
    let owners = owners(project, &owned);
    let mut diagnostics = Vec::new();
    cleanups(project, &owned, &released, &mut diagnostics);
    for procedure in procedures {
        procedure.diagnostics(&owners, &mut diagnostics);
    }
    diagnostics
}
