use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    sync::LazyLock,
};

use regex::Regex;
use tree_sitter::Node;

use crate::{
    config::Config,
    diagnostic::Diagnostic,
    flow::{self, Edge, Flow, Step, unwrap},
    library_calls::{
        ArgumentSelector, Call, Family, Part, RequestAction, Storage, arguments, assigned,
        body_nodes, call_name, declarator, escapes, executable, identity, int, operator, read_call,
        status_tests, storage,
    },
    memory::formal_arguments,
    semantic::Project,
    source::{Unit, child, children, fields, walk},
};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Unchecked {
    site: usize,
    carriers: BTreeSet<Storage>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct RequestOperation {
    action: RequestAction,
    place: Storage,
    array: bool,
    whole_array: bool,
    count: Option<i64>,
    collective: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Request {
    site: usize,
    operation: RequestOperation,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Completion {
    flag: Storage,
    request: Request,
    site: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Split {
    site: usize,
    begin: String,
    arguments: Vec<(Option<String>, Option<String>)>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct State {
    values: BTreeSet<(Storage, usize)>,
    unchecked: BTreeSet<Unchecked>,
    ranks: BTreeSet<(Storage, String)>,
    active: BTreeSet<Request>,
    persistent: BTreeSet<Request>,
    confirmed_active: BTreeSet<Request>,
    confirmed_persistent: BTreeSet<Request>,
    tested: BTreeSet<Completion>,
    splits: BTreeSet<Split>,
}

impl State {
    fn join(left: &Self, right: &Self) -> Self {
        Self {
            values: left.values.union(&right.values).cloned().collect(),
            unchecked: left.unchecked.union(&right.unchecked).cloned().collect(),
            ranks: left.ranks.union(&right.ranks).cloned().collect(),
            active: left.active.union(&right.active).cloned().collect(),
            persistent: left.persistent.union(&right.persistent).cloned().collect(),
            confirmed_active: left
                .confirmed_active
                .intersection(&right.confirmed_active)
                .cloned()
                .collect(),
            confirmed_persistent: left
                .confirmed_persistent
                .intersection(&right.confirmed_persistent)
                .cloned()
                .collect(),
            tested: left.tested.intersection(&right.tested).cloned().collect(),
            splits: left.splits.union(&right.splits).cloned().collect(),
        }
    }

    fn retain_confirmation(&mut self) {
        self.confirmed_active
            .retain(|item| self.active.contains(item));
        self.confirmed_persistent
            .retain(|item| self.persistent.contains(item));
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CollectiveEvent {
    routine: String,
    parameters: Vec<Option<String>>,
}

type StatusValues = BTreeSet<(Storage, usize)>;
type StatusEffect = (StatusValues, BTreeSet<Unchecked>, BTreeSet<usize>);

struct Procedure<'a, 'p> {
    project: &'p Project<'a>,
    config: &'p Config,
    scope: usize,
    unit: &'p Unit,
    calls: BTreeMap<usize, Call<'a>>,
    origins: BTreeMap<usize, Node<'a>>,
    requests: BTreeMap<usize, (RequestOperation, Option<Storage>)>,
    statuses: BTreeMap<usize, (bool, Option<Node<'a>>)>,
    writes: BTreeMap<usize, BTreeSet<Storage>>,
    collectives: BTreeMap<usize, (String, CollectiveEvent)>,
    splits: BTreeMap<usize, (Split, bool)>,
}

impl<'a, 'p> Procedure<'a, 'p> {
    fn new(project: &'p Project<'a>, config: &'p Config, scope: usize, nodes: &[Node<'a>]) -> Self {
        let unit = &project.units[project.scopes[scope].unit];
        let mut procedure = Self {
            project,
            config,
            scope,
            unit,
            calls: BTreeMap::new(),
            origins: BTreeMap::new(),
            requests: BTreeMap::new(),
            statuses: BTreeMap::new(),
            writes: BTreeMap::new(),
            collectives: BTreeMap::new(),
            splits: BTreeMap::new(),
        };
        for node in nodes
            .iter()
            .copied()
            .filter(|node| node.kind() == "subroutine_call")
        {
            let site = node.id();
            let Some(call) = read_call(project, config, scope, node) else {
                continue;
            };
            let selector = config.safety.expected_status(&call);
            let status = selector.and_then(|selector| call.argument(selector));
            procedure
                .statuses
                .insert(site, (selector.is_some(), status));
            let mut writes = BTreeSet::new();
            if let Some(status) = status {
                writes.extend(procedure.place(status));
            }
            if let Some(signature) = config.safety.mpi(&call) {
                writes.extend(
                    signature
                        .outputs
                        .iter()
                        .filter_map(|arg| procedure.actual(&call, arg)),
                );
                if let Some(rank) = signature.rank.as_ref().filter(|rank| rank.enabled) {
                    writes.extend(procedure.actual(&call, &rank.rank));
                }
                if let Some(collective) = signature
                    .collective
                    .as_ref()
                    .filter(|collective| collective.enabled)
                    && let Some(comm) = call.argument(&collective.communicator)
                {
                    let event = CollectiveEvent {
                        routine: call.name.clone(),
                        parameters: collective
                            .compare_arguments
                            .iter()
                            .map(|selector| {
                                call.argument(selector)
                                    .and_then(|node| procedure.comparable(node))
                            })
                            .collect(),
                    };
                    procedure
                        .collectives
                        .insert(site, (procedure.identity(comm), event));
                }
                if let Some(request) = signature.request.as_ref().filter(|request| request.enabled)
                    && let Some(place) = procedure.actual(&call, &request.handle)
                {
                    let count = request
                        .count
                        .as_ref()
                        .and_then(|selector| call.argument(selector));
                    let whole_array =
                        count.is_some_and(|count| procedure.is_size_of(count, &place));
                    let count = count.and_then(|node| int(procedure.unit.text(node)));
                    let flag = request
                        .completion_flag
                        .as_ref()
                        .and_then(|selector| procedure.actual(&call, selector));
                    writes.extend(flag.clone());
                    if matches!(
                        request.action,
                        RequestAction::Post
                            | RequestAction::Initialize
                            | RequestAction::Start
                            | RequestAction::Wait
                            | RequestAction::Test
                            | RequestAction::Free
                    ) {
                        writes.insert(place.clone());
                    }
                    procedure.requests.insert(
                        site,
                        (
                            RequestOperation {
                                action: request.action,
                                place,
                                array: request.count.is_some(),
                                whole_array,
                                count,
                                collective: request.collective,
                            },
                            flag,
                        ),
                    );
                }
            }
            if let Some(split) = procedure.split(&call, site, status) {
                procedure.splits.insert(site, split);
            }
            procedure.writes.insert(site, writes);
            procedure.origins.insert(site, node);
            procedure.calls.insert(site, call);
        }
        procedure
    }

    fn place(&self, node: Node<'_>) -> Option<Storage> {
        storage(self.project, self.scope, node)
    }

    fn identity(&self, node: Node<'_>) -> String {
        identity(self.project, self.scope, node)
    }

    fn actual(&self, call: &Call<'a>, argument: &ArgumentSelector) -> Option<Storage> {
        call.argument(argument).and_then(|node| self.place(node))
    }

    fn comparable(&self, node: Node<'a>) -> Option<String> {
        let node = unwrap(node);
        let name = self.identity(node);
        if int(&name).is_some()
            || name.starts_with("mpi_")
            || self
                .project
                .expression(self.project.scopes[self.scope].unit, node)
                .is_some_and(|id| self.project.symbols[id].attrs.contains("parameter"))
        {
            return Some(name);
        }
        None
    }

    fn is_size_of(&self, node: Node<'a>, place: &Storage) -> bool {
        let node = unwrap(node);
        if node.kind() != "call_expression"
            || node
                .named_child(0)
                .is_none_or(|node| self.unit.name(node) != "size")
        {
            return false;
        }
        let args = arguments(self.unit, node);
        args.len() == 1 && self.place(args[0].node).as_ref() == Some(place)
    }

    fn split(
        &self,
        call: &Call<'a>,
        site: usize,
        status: Option<Node<'a>>,
    ) -> Option<(Split, bool)> {
        if call.family != Family::Petsc {
            return None;
        }
        let policy = &self.config.safety.petsc_splits;
        let (begin, pair, start) =
            if let Some(pair) = policy.get(&call.name).filter(|pair| pair.enabled) {
                (&call.name, pair, true)
            } else {
                let (begin, pair) = policy
                    .iter()
                    .find(|(_, pair)| pair.enabled && pair.end == call.name)?;
                (begin, pair, false)
            };
        let mut arguments: Vec<_> = if pair.context.is_empty() {
            call.arguments
                .iter()
                .filter(|arg| Some(arg.node) != status)
                .map(|arg| (arg.keyword.clone(), Some(self.identity(arg.node))))
                .collect()
        } else {
            pair.context
                .iter()
                .map(|selector| {
                    (
                        None,
                        call.argument(selector).map(|node| self.identity(node)),
                    )
                })
                .collect()
        };
        if pair.context.is_empty() {
            arguments.sort_by(|left, right| left.0.cmp(&right.0));
        }
        Some((
            Split {
                site,
                begin: begin.clone(),
                arguments,
            },
            start,
        ))
    }

    fn diagnostic(
        &self,
        code: &str,
        message: impl Into<String>,
        node: Node<'_>,
        certain: bool,
    ) -> Diagnostic {
        let diagnostic = Diagnostic::new(code, message, self.unit.location(node))
            .with_subject(&self.project.scopes[self.scope].qualified);
        if certain {
            diagnostic
        } else {
            diagnostic.warning()
        }
    }

    fn status_diagnostic(&self, site: usize, message: &str, certain: bool) -> Diagnostic {
        let call = &self.calls[&site];
        let code = if call.family == Family::Mpi {
            "PRL001"
        } else {
            "PRL002"
        };
        self.diagnostic(
            code,
            format!("{} {message}", call.name.to_ascii_uppercase()),
            self.origins[&site],
            certain,
        )
    }

    fn references(&self, node: Node<'_>) -> BTreeSet<Storage> {
        walk(node).filter_map(|node| self.place(node)).collect()
    }

    fn observed(&self, node: Node<'a>) -> BTreeSet<Storage> {
        if matches!(
            node.kind(),
            "selector" | "parenthesized_expression" | "write_statement" | "print_statement"
        ) {
            return self.references(node);
        }
        if let Some((_, right)) = assigned(node) {
            let value = unwrap(right);
            if matches!(value.kind(), "logical_expression" | "relational_expression") {
                return self.references(right);
            }
        }
        if node.kind() != "subroutine_call" {
            return BTreeSet::new();
        }
        if let Some(call) = self.calls.get(&node.id()) {
            return self
                .config
                .safety
                .status_handlers
                .get(&call.name)
                .and_then(|selector| self.actual(call, selector))
                .into_iter()
                .collect();
        }
        if let Some(name) = call_name(self.project, self.scope, node)
            && let Some(selector) = self
                .config
                .safety
                .status_handlers
                .get(&name)
                .filter(|_| name != "mpi_abort" && name != "mpiu_abort")
        {
            let call = Call {
                family: Family::Petsc,
                name,
                arguments: arguments(self.unit, node),
            };
            return self.actual(&call, selector).into_iter().collect();
        }
        formal_arguments(self.project, self.scope, node)
            .into_iter()
            .filter_map(|(actual, formal)| {
                self.project.symbols[formal]
                    .attrs
                    .contains("intent(in)")
                    .then(|| self.place(actual))?
            })
            .collect()
    }

    fn writes(&self, node: Node<'a>, destination: Option<&Storage>) -> BTreeSet<Storage> {
        let mut writes = self.writes.get(&node.id()).cloned().unwrap_or_default();
        writes.extend(destination.cloned());
        if node.kind() == "subroutine_call" && !self.calls.contains_key(&node.id()) {
            writes.extend(
                formal_arguments(self.project, self.scope, node)
                    .into_iter()
                    .filter_map(|(actual, formal)| {
                        (!self.project.symbols[formal].attrs.contains("intent(in)"))
                            .then(|| self.place(actual))?
                    }),
            );
        }
        writes
    }

    fn status_effect(
        &self,
        node: Node<'a>,
        facts: &State,
        destination: Option<&Storage>,
        source: Option<&Storage>,
        writes: &BTreeSet<Storage>,
    ) -> StatusEffect {
        let observed = self.observed(node);
        let mut values: BTreeSet<_> = facts
            .values
            .iter()
            .filter(|(place, _)| !written(writes, place))
            .cloned()
            .collect();
        if let (Some(destination), Some(source)) = (destination, source) {
            values.extend(
                facts
                    .values
                    .iter()
                    .filter(|(place, _)| place == source)
                    .map(|(_, site)| (destination.clone(), *site)),
            );
        }
        let mut unchecked = BTreeSet::new();
        let mut lost = BTreeSet::new();
        for status in &facts.unchecked {
            if !status.carriers.is_disjoint(&observed) {
                continue;
            }
            let mut carriers: BTreeSet<_> = status
                .carriers
                .iter()
                .filter(|place| !written(writes, place))
                .cloned()
                .collect();
            if let (Some(destination), Some(source)) = (destination, source)
                && status.carriers.contains(source)
            {
                carriers.insert(destination.clone());
            }
            if carriers.is_empty() {
                lost.insert(status.site);
            } else {
                unchecked.insert(Unchecked {
                    site: status.site,
                    carriers,
                });
            }
        }
        if !self.config.safety.terminates(self.calls.get(&node.id()))
            && let Some(output) = self
                .statuses
                .get(&node.id())
                .and_then(|(_, node)| *node)
                .and_then(|node| self.place(node))
        {
            values.insert((output.clone(), node.id()));
            unchecked.insert(Unchecked {
                site: node.id(),
                carriers: BTreeSet::from([output]),
            });
        }
        (values, unchecked, lost)
    }

    fn transfer(&self, step: &Step<'a>, before: &State) -> State {
        let Some(node) = step.node else {
            return before.clone();
        };
        let site = node.id();
        let assignment = assigned(node).filter(|_| node.kind() == "assignment_statement");
        let destination = assignment.and_then(|(left, _)| self.place(left));
        let source = assignment.and_then(|(_, right)| self.place(right));
        let writes = self.writes(node, destination.as_ref());
        let (values, unchecked, _) =
            self.status_effect(node, before, destination.as_ref(), source.as_ref(), &writes);
        let mut after = before.clone();
        after.values = values;
        after.unchecked = unchecked;
        after
            .ranks
            .retain(|(storage, _)| !written(&writes, storage));
        after.tested.retain(|item| !written(&writes, &item.flag));
        if let Some(destination) = &destination {
            if !before.ranks.is_empty()
                && let Some((_, right)) = assignment
            {
                let dependencies = self.references(right);
                after.ranks.extend(
                    before
                        .ranks
                        .iter()
                        .filter(|(place, _)| dependencies.contains(place))
                        .map(|(_, communicator)| (destination.clone(), communicator.clone())),
                );
            }
            if let Some(source) = &source {
                after.tested.extend(
                    before
                        .tested
                        .iter()
                        .filter(|item| &item.flag == source)
                        .map(|item| Completion {
                            flag: destination.clone(),
                            ..item.clone()
                        }),
                );
            }
        }
        for (_, replaced) in self.request_replacements(node) {
            after
                .active
                .retain(|item| !contains(&replaced, &item.operation.place));
            after
                .persistent
                .retain(|item| !contains(&replaced, &item.operation.place));
        }
        if let Some(call) = self.calls.get(&site)
            && let Some(signature) = self
                .config
                .safety
                .mpi(call)
                .and_then(|signature| signature.rank.as_ref())
                .filter(|rank| rank.enabled)
            && let (Some(rank), Some(comm)) = (
                self.actual(call, &signature.rank),
                call.argument(&signature.communicator),
            )
        {
            after.ranks.insert((rank, self.identity(comm)));
        }
        if let Some((operation, flag)) = self.requests.get(&site) {
            self.request_effect(site, operation, flag.as_ref(), &mut after);
        }
        if let Some((split, start)) = self.splits.get(&site) {
            if *start {
                after.splits.insert(split.clone());
            } else {
                after
                    .splits
                    .retain(|item| item.begin != split.begin || item.arguments != split.arguments);
            }
        }
        after.retain_confirmation();
        after
    }

    fn request_replacements(&self, node: Node<'a>) -> Vec<(Node<'a>, Storage)> {
        if node.kind() == "assignment_statement" {
            if let Some((left, right)) = assigned(node)
                && let Some(place) = self
                    .place(left)
                    .filter(|place| Some(place) != self.place(right).as_ref())
            {
                return vec![(left, place)];
            }
        } else if node.kind() == "deallocate_statement" {
            return children(node)
                .into_iter()
                .filter_map(|target| Some((target, self.place(target)?)))
                .collect();
        }
        Vec::new()
    }

    fn request_effect(
        &self,
        site: usize,
        operation: &RequestOperation,
        flag: Option<&Storage>,
        state: &mut State,
    ) {
        let matches = |item: &Request| covers(self.project, operation, &item.operation);
        let created = Request {
            site,
            operation: operation.clone(),
        };
        let certain = self
            .statuses
            .get(&site)
            .is_none_or(|(_, status)| status.is_none());
        match operation.action {
            RequestAction::Post => {
                state.active.retain(|item| !matches(item));
                state.tested.retain(|item| !matches(&item.request));
                if certain {
                    state.confirmed_active.insert(created.clone());
                }
                state.active.insert(created);
            }
            RequestAction::Initialize => {
                state.persistent.retain(|item| !matches(item));
                if certain {
                    state.confirmed_persistent.insert(created.clone());
                }
                state.persistent.insert(created);
            }
            RequestAction::Start => {
                let mut sources: Vec<_> = state
                    .persistent
                    .iter()
                    .filter(|item| matches(item))
                    .map(|item| Request {
                        site,
                        operation: item.operation.clone(),
                    })
                    .collect();
                if sources.is_empty() || operation.array {
                    sources.push(created);
                }
                state.tested.retain(|item| !matches(&item.request));
                if certain {
                    state.confirmed_active.extend(sources.iter().cloned());
                }
                state.active.extend(sources);
            }
            RequestAction::Wait | RequestAction::Free => {
                state.active.retain(|item| !matches(item));
                state.tested.retain(|item| !matches(&item.request));
                if operation.action == RequestAction::Free {
                    state.persistent.retain(|item| !matches(item));
                }
            }
            RequestAction::Test => {
                if let Some(flag) = flag {
                    state.tested.retain(|item| &item.flag != flag);
                    state
                        .tested
                        .extend(
                            state
                                .active
                                .iter()
                                .filter(|item| matches(item))
                                .map(|item| Completion {
                                    flag: flag.clone(),
                                    request: item.clone(),
                                    site,
                                }),
                        );
                }
            }
            RequestAction::Cancel => {}
        }
    }

    fn refine(&self, edge: &Edge<'a>, before: &State) -> Option<State> {
        let mut state = before.clone();
        for (status, success) in status_tests(self.project, self.config, self.scope, edge) {
            let origins: BTreeSet<_> = state
                .values
                .iter()
                .filter(|(place, _)| *place == status)
                .map(|(_, site)| *site)
                .collect();
            if origins.len() != 1 {
                continue;
            }
            if success {
                state.confirmed_active.extend(
                    state
                        .active
                        .iter()
                        .filter(|item| origins.contains(&item.site))
                        .cloned(),
                );
                state.confirmed_persistent.extend(
                    state
                        .persistent
                        .iter()
                        .filter(|item| origins.contains(&item.site))
                        .cloned(),
                );
            } else {
                state.active.retain(|item| !origins.contains(&item.site));
                state
                    .persistent
                    .retain(|item| !origins.contains(&item.site));
                state.tested.retain(|item| !origins.contains(&item.site));
                state.splits.retain(|item| !origins.contains(&item.site));
            }
        }
        if let Some(condition) = edge.condition {
            let mut completed_flags = BTreeSet::new();
            true_flags(self, condition, edge.truth, &mut completed_flags);
            let completed: BTreeSet<_> = state
                .tested
                .iter()
                .filter(|item| completed_flags.contains(&item.flag))
                .map(|item| &item.request)
                .collect();
            state.active.retain(|item| !completed.contains(item));
        }
        state.retain_confirmation();
        Some(state)
    }

    fn unfinished(&self, facts: &State, finalization: bool, diagnostics: &mut Vec<Diagnostic>) {
        let end = if finalization {
            "library finalization"
        } else {
            "normal exit"
        };
        for request in facts.active.union(&facts.persistent) {
            let place = &request.operation.place;
            if !finalization && escapes(self.project, self.scope, place) {
                continue;
            }
            let symbol = &self.project.symbols[place.root];
            let message = if facts.active.contains(request) {
                format!("Request '{}' may remain active at {end}", symbol.name)
            } else {
                format!(
                    "Persistent request '{}' may remain unfreed at {end}",
                    symbol.name
                )
            };
            diagnostics.push(self.diagnostic(
                "PRL003",
                message,
                self.origins[&request.site],
                false,
            ));
        }
        for split in &facts.splits {
            let routine = split.begin.to_ascii_uppercase();
            diagnostics.push(self.diagnostic(
                "PRL003",
                format!("{routine} has no matching end call before {end}"),
                self.origins[&split.site],
                false,
            ));
        }
    }

    fn status_outputs(&self, place: &Storage) -> bool {
        let symbol = &self.project.symbols[place.root];
        let info = &self.project.scopes[self.scope];
        symbol.scope == self.scope
            && (info.result.as_ref() == Some(&symbol.name)
                || (info.parameters.contains(&symbol.name)
                    && (symbol.attrs.contains("intent(inout)")
                        || symbol.attrs.contains("intent(out)"))))
    }

    fn request_checks(
        &self,
        node: Node<'a>,
        before: &State,
        operation: &RequestOperation,
        diagnostics: &mut Vec<Diagnostic>,
    ) {
        let overlaps = |item: &&Request| overlaps(self.project, operation, &item.operation);
        let active: Vec<_> = before.active.iter().filter(overlaps).collect();
        let initialized: Vec<_> = before.persistent.iter().filter(overlaps).collect();
        let confirmed = |item: &&Request| {
            before.confirmed_active.contains(*item) || before.confirmed_persistent.contains(*item)
        };
        let request = &self.project.symbols[operation.place.root].name;
        let (message, certain) = match operation.action {
            RequestAction::Post | RequestAction::Start | RequestAction::Initialize
                if !active.is_empty() =>
            {
                (
                    format!("Operation may reuse active request '{request}'"),
                    active.iter().any(confirmed),
                )
            }
            RequestAction::Post | RequestAction::Initialize if !initialized.is_empty() => (
                format!("Operation can overwrite initialized persistent request '{request}'"),
                initialized.iter().any(confirmed),
            ),
            RequestAction::Cancel
                if active
                    .iter()
                    .chain(&initialized)
                    .any(|item| item.operation.collective) =>
            {
                (
                    "MPI_CANCEL is invalid for a collective request".into(),
                    active
                        .iter()
                        .chain(&initialized)
                        .filter(|item| item.operation.collective)
                        .any(confirmed),
                )
            }
            RequestAction::Free if !active.is_empty() => (
                if active.iter().any(|item| item.operation.collective) {
                    "MPI_REQUEST_FREE is invalid for an active collective request".into()
                } else {
                    format!("Freeing active request '{request}' can hide unfinished communication")
                },
                active
                    .iter()
                    .filter(|item| item.operation.collective)
                    .any(confirmed),
            ),
            _ => return,
        };
        diagnostics.push(self.diagnostic("PRL004", message, node, certain));
    }

    fn diagnostics(&self, control: &Flow<'a>) -> Vec<Diagnostic> {
        let facts = flow::solve(
            control,
            State::default(),
            |step, facts| self.transfer(step, facts),
            State::join,
            |edge, facts| self.refine(edge, facts),
        );
        let mut result = Vec::new();
        for (index, (step, before)) in control.steps.iter().zip(&facts).enumerate() {
            let (Some(node), Some(before)) = (step.node, before) else {
                continue;
            };
            let site = node.id();
            if matches!(self.statuses.get(&site), Some((true, None))) {
                result.push(self.status_diagnostic(
                    site,
                    "requires an explicit error-status argument",
                    true,
                ));
            }
            let assignment = assigned(node).filter(|_| node.kind() == "assignment_statement");
            let destination = assignment.and_then(|(left, _)| self.place(left));
            let source = assignment.and_then(|(_, right)| self.place(right));
            let writes = self.writes(node, destination.as_ref());
            let (_, unchecked, lost) =
                self.status_effect(node, before, destination.as_ref(), source.as_ref(), &writes);
            for dropped in lost {
                result.push(self.status_diagnostic(
                    dropped,
                    "error status may be overwritten before it is checked",
                    false,
                ));
            }
            if step.edges.is_empty() {
                for item in unchecked {
                    if !item.carriers.iter().any(|place| self.status_outputs(place)) {
                        result.push(self.status_diagnostic(
                            item.site,
                            "error status may remain unchecked at termination",
                            false,
                        ));
                    }
                }
            }
            if let Some((operation, _)) = self.requests.get(&site) {
                self.request_checks(node, before, operation, &mut result);
            }
            for (target, destination) in self.request_replacements(node) {
                let replaced = before
                    .active
                    .union(&before.persistent)
                    .any(|item| contains(&destination, &item.operation.place));
                if replaced {
                    let certain = before
                        .confirmed_active
                        .union(&before.confirmed_persistent)
                        .any(|item| contains(&destination, &item.operation.place));
                    let symbol = &self.project.symbols[destination.root];
                    let null = assignment.is_some_and(|(_, rhs)| {
                        self.config
                            .safety
                            .request_nulls
                            .contains(&self.identity(rhs))
                    });
                    let action = if null {
                        "Resetting"
                    } else if node.kind() == "deallocate_statement" {
                        "Deallocating"
                    } else {
                        "Overwriting"
                    };
                    result.push(self.diagnostic(
                        "PRL004",
                        format!(
                            "{action} request '{}' can discard unfinished communication",
                            symbol.name
                        ),
                        target,
                        certain,
                    ));
                }
            }
            if self
                .calls
                .get(&site)
                .is_some_and(|call| self.config.safety.finalizers.contains(&call.name))
            {
                self.unfinished(before, true, &mut result);
            }
            if step.edges.len() > 1
                && !before.ranks.is_empty()
                && self.divergent(control, index, node, before)
            {
                result.push(self.diagnostic("PRL005", "MPI ranks may use different collective calls or arguments on the same communicator",
                    node, false));
            }
        }
        if let Some(exit) = &facts[control.exit] {
            for status in &exit.unchecked {
                if !status
                    .carriers
                    .iter()
                    .any(|place| self.status_outputs(place))
                {
                    result.push(self.status_diagnostic(
                        status.site,
                        "error status may remain unchecked at normal exit",
                        false,
                    ));
                }
            }
            self.unfinished(exit, false, &mut result);
        }
        result
    }

    fn divergent(&self, flow: &Flow<'a>, branch: usize, node: Node<'a>, facts: &State) -> bool {
        let dependencies = self.references(node);
        let candidates: BTreeSet<_> = facts
            .ranks
            .iter()
            .filter(|(place, comm)| {
                dependencies.contains(place)
                    && !self.config.safety.self_communicators.contains(comm)
            })
            .map(|(_, comm)| comm)
            .collect();
        if candidates.is_empty() {
            return false;
        }
        let targets: Vec<_> = flow.steps[branch]
            .edges
            .iter()
            .map(|edge| edge.target)
            .collect();
        let join = rejoin(flow, &targets, branch);
        for comm in candidates {
            let events: BTreeMap<_, _> = self
                .collectives
                .iter()
                .filter(|(_, (id, _))| id == comm)
                .map(|(site, (_, event))| (*site, event))
                .collect();
            if events.is_empty() {
                continue;
            }
            let comm_selector = ArgumentSelector::new(0, &["comm"]);
            let aborts: BTreeSet<_> = self
                .calls
                .iter()
                .filter(|(_, call)| {
                    self.config.safety.terminates(Some(call))
                        && call
                            .argument(&comm_selector)
                            .is_some_and(|node| self.identity(node) == *comm)
                })
                .map(|(site, _)| *site)
                .collect();
            let traces: Vec<_> = targets
                .iter()
                .map(|target| {
                    traces(
                        flow,
                        *target,
                        join,
                        branch,
                        &events,
                        &aborts,
                        self.config.safety.collective_path_limit,
                    )
                })
                .collect();
            if traces.iter().all(|path| !path.is_empty())
                && traces
                    .iter()
                    .skip(1)
                    .any(|path| !same_traces(&traces[0], path))
            {
                return true;
            }
        }
        false
    }
}

fn true_flags(
    procedure: &Procedure<'_, '_>,
    node: Node<'_>,
    truth: bool,
    result: &mut BTreeSet<Storage>,
) {
    let node = unwrap(node);
    if procedure.unit.name(node).starts_with(".not.") {
        if let Some(argument) = node
            .child_by_field_name("argument")
            .or_else(|| node.named_child(0))
        {
            true_flags(procedure, argument, !truth, result);
        }
        return;
    }
    if let (Some(left), Some(right)) = (
        node.child_by_field_name("left"),
        node.child_by_field_name("right"),
    ) {
        let op = operator(procedure.unit, left, right);
        if (op == ".and." && truth) || (op == ".or." && !truth) {
            true_flags(procedure, left, truth, result);
            true_flags(procedure, right, truth, result);
            return;
        }
        if matches!(op.as_str(), ".eqv." | ".neqv.") {
            for (value, literal) in [(left, right), (right, left)] {
                if literal.kind() == "boolean_literal" {
                    let equal = truth == (op == ".eqv.");
                    true_flags(
                        procedure,
                        value,
                        equal == procedure.unit.name(literal).starts_with(".true."),
                        result,
                    );
                    return;
                }
            }
        }
    }
    if truth {
        result.extend(procedure.place(node));
    }
}

fn contains(outer: &Storage, inner: &Storage) -> bool {
    outer.root == inner.root && inner.parts.starts_with(&outer.parts)
}

fn written(writes: &BTreeSet<Storage>, place: &Storage) -> bool {
    writes.iter().any(|target| contains(target, place))
}

fn covers(project: &Project<'_>, outer: &RequestOperation, inner: &RequestOperation) -> bool {
    if outer.place == inner.place {
        if !inner.array {
            return !outer.array || outer.whole_array || outer.count.is_some_and(|count| count > 0);
        }
        return outer.array
            && (outer.whole_array
                || outer
                    .count
                    .zip(inner.count)
                    .is_some_and(|(outer, inner)| outer >= inner));
    }
    if !outer.array
        || !contains(&outer.place, &inner.place)
        || inner.place.parts.len() != outer.place.parts.len() + 1
    {
        return false;
    }
    let Some(Part::Indices(indices)) = inner.place.parts.last() else {
        return false;
    };
    if outer.whole_array {
        return true;
    }
    let [Some(index)] = indices.as_slice() else {
        return false;
    };
    if !outer.place.parts.is_empty() {
        return false;
    }
    let symbol = &project.symbols[outer.place.root];
    let unit = &project.units[project.scopes[symbol.scope].unit];
    let extent = walk(declarator(symbol))
        .find(|node| node.kind() == "size")
        .and_then(|node| node.named_child(0));
    let lower = extent.map(|node| unit.name(node)).and_then(|text| {
        text.split_once(':').map_or(
            Some(1),
            |(lo, hi)| if hi.is_empty() { None } else { int(lo) },
        )
    });
    outer.count.zip(lower).is_some_and(|(count, lower)| {
        *index >= lower && index.checked_sub(lower).is_some_and(|d| d < count)
    })
}

fn overlaps(project: &Project<'_>, left: &RequestOperation, right: &RequestOperation) -> bool {
    covers(project, left, right) || covers(project, right, left)
}

fn rejoin(flow: &Flow<'_>, targets: &[usize], branch: usize) -> Option<usize> {
    let distances: Vec<_> = targets
        .iter()
        .map(|start| {
            let mut distances = BTreeMap::new();
            let mut pending = VecDeque::from([(*start, 0usize)]);
            while let Some((point, depth)) = pending.pop_front() {
                if point == branch || distances.contains_key(&point) {
                    continue;
                }
                distances.insert(point, depth);
                pending.extend(
                    flow.steps[point]
                        .edges
                        .iter()
                        .map(|edge| (edge.target, depth + 1)),
                );
            }
            distances
        })
        .collect();
    distances
        .first()?
        .keys()
        .filter_map(|point| {
            let depths: Option<Vec<_>> = distances
                .iter()
                .map(|lengths| lengths.get(point).copied())
                .collect();
            let depths = depths?;
            Some((*depths.iter().max()?, depths.iter().sum::<usize>(), *point))
        })
        .min()
        .map(|(_, _, point)| point)
}

fn traces<'a>(
    flow: &Flow<'_>,
    start: usize,
    join: Option<usize>,
    branch: usize,
    calls: &BTreeMap<usize, &'a CollectiveEvent>,
    aborts: &BTreeSet<usize>,
    limit: usize,
) -> BTreeSet<Vec<&'a CollectiveEvent>> {
    let mut result = BTreeSet::new();
    let mut pending = vec![(start, Vec::new(), BTreeSet::new())];
    for _ in 0..limit {
        let Some((point, mut path, mut visited)) = pending.pop() else {
            return result;
        };
        if Some(point) == join || point == branch || !visited.insert(point) {
            result.insert(path);
            continue;
        }
        let step = &flow.steps[point];
        let site = step.node.map(|node| node.id());
        if let Some(event) = site.and_then(|site| calls.get(&site)) {
            path.push(*event);
        }
        if !step.edges.is_empty() {
            pending.extend(
                step.edges
                    .iter()
                    .map(|edge| (edge.target, path.clone(), visited.clone())),
            );
        } else if site.is_none_or(|site| !aborts.contains(&site)) {
            result.insert(path);
        }
    }
    if pending.is_empty() {
        result
    } else {
        BTreeSet::new()
    }
}

fn same_traces(
    left: &BTreeSet<Vec<&CollectiveEvent>>,
    right: &BTreeSet<Vec<&CollectiveEvent>>,
) -> bool {
    let matches = |first: &[&CollectiveEvent], second: &[&CollectiveEvent]| {
        first.len() == second.len()
            && first.iter().zip(second).all(|(left, right)| {
                left.routine == right.routine
                    && left
                        .parameters
                        .iter()
                        .zip(&right.parameters)
                        .all(|(first, second)| {
                            first.is_none() || second.is_none() || first == second
                        })
            })
    };
    left.iter()
        .all(|path| right.iter().any(|other| matches(path, other)))
        && right
            .iter()
            .all(|path| left.iter().any(|other| matches(path, other)))
}

#[derive(Default)]
struct Region {
    first: usize,
    last: usize,
    private: BTreeSet<String>,
    indices: BTreeSet<String>,
    protect: bool,
}

impl Region {
    fn contains(&self, node: Node<'_>) -> bool {
        (self.first..self.last).contains(&node.start_byte())
    }
}

fn clause_names(clauses: &str, pattern: &Regex) -> BTreeSet<String> {
    pattern
        .captures_iter(clauses)
        .flat_map(|capture| {
            capture
                .get(1)
                .map(|part| {
                    let names = part.as_str().rsplit(':').next().unwrap_or_default();
                    names.split(',').map(str::trim)
                })
                .into_iter()
                .flatten()
        })
        .map(str::to_owned)
        .collect()
}

fn next_statement(node: Node<'_>) -> Option<Node<'_>> {
    let mut sibling = node.next_named_sibling();
    while let Some(node) = sibling {
        if node.kind() != "comment" {
            return Some(node);
        }
        sibling = node.next_named_sibling();
    }
    None
}

fn openmp_regions(unit: &Unit, nodes: &[Node<'_>]) -> (Vec<Region>, Vec<Region>) {
    static CLAUSES: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"\b(?:private|firstprivate|lastprivate|reduction|linear)\s*\(([^)]*)\)")
            .unwrap()
    });
    static ALL_PRIVATE: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"\bdefault\s*\(\s*(?:private|firstprivate)\s*\)").unwrap());
    let mut directives: Vec<(Node<'_>, String)> = Vec::new();
    for node in nodes
        .iter()
        .copied()
        .filter(|node| node.kind() == "comment")
    {
        let raw = unit.text(node).trim_start().to_ascii_lowercase();
        let Some(text) = raw.strip_prefix("!$omp") else {
            continue;
        };
        let text = text
            .trim_start()
            .trim_start_matches('&')
            .split('!')
            .next()
            .unwrap_or_default()
            .trim();
        if let Some((_, previous)) = directives
            .last_mut()
            .filter(|(_, text)| text.ends_with('&'))
        {
            previous.pop();
            previous.push(' ');
            previous.push_str(text);
        } else {
            directives.push((node, text.to_owned()));
        }
    }
    let mut active = Vec::new();
    let mut protected = Vec::new();
    let mut open: Vec<(&str, Region)> = Vec::new();
    for (node, directive) in &directives {
        let mut tokens = directive
            .split(|c: char| c.is_whitespace() || c == '(')
            .filter(|item| !item.is_empty());
        let Some(first) = tokens.next() else {
            continue;
        };
        match first {
            "end" => {
                let Some(last) = tokens.next() else {
                    continue;
                };
                if let Some(index) = open.iter().rposition(|(name, _)| *name == last) {
                    let (name, mut region) = open.remove(index);
                    region.last = node.start_byte();
                    if name == "parallel" {
                        active.push(region);
                    } else {
                        protected.push(region);
                    }
                }
            }
            "parallel" | "do" | "loop" => {
                let mut private = clause_names(directive, &CLAUSES);
                if ALL_PRIVATE.is_match(directive) {
                    private.insert("*".into());
                }
                if first != "parallel" || matches!(tokens.next(), Some("do" | "loop")) {
                    if let Some(body) =
                        next_statement(*node).filter(|body| body.kind() == "do_loop")
                    {
                        let region = Region {
                            first: body.start_byte(),
                            last: body.end_byte(),
                            private,
                            ..Region::default()
                        };
                        if first == "parallel" {
                            active.push(region);
                        } else {
                            protected.push(region);
                        }
                    }
                } else {
                    open.push((
                        first,
                        Region {
                            first: node.end_byte(),
                            private,
                            ..Region::default()
                        },
                    ));
                }
            }
            "critical" | "ordered" | "single" | "master" | "masked" => {
                open.push((
                    first,
                    Region {
                        first: node.end_byte(),
                        protect: true,
                        ..Region::default()
                    },
                ));
            }
            "atomic" => {
                if let Some(body) = next_statement(*node) {
                    protected.push(Region {
                        first: body.start_byte(),
                        last: body.end_byte(),
                        protect: true,
                        ..Region::default()
                    });
                }
            }
            _ => {}
        }
    }
    (active, protected)
}

fn concurrent_regions(unit: &Unit, nodes: &[Node<'_>]) -> Vec<Region> {
    static LOCAL: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"\b(?:local|local_init|reduce)\s*\(([^)]*)\)").unwrap());
    static CONTROL: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"(?:\(|,)\s*([a-z][a-z0-9_]*)\s*=").unwrap());
    nodes
        .iter()
        .copied()
        .filter_map(|node| {
            if node.kind() != "do_loop" {
                return None;
            }
            let header = child(node, "do_statement")?;
            let concurrent = child(header, "concurrent_statement")?;
            let raw = unit.text(concurrent).to_ascii_lowercase();
            let private = clause_names(&raw, &LOCAL);
            let indices = CONTROL
                .captures_iter(&raw)
                .map(|capture| capture[1].to_owned())
                .collect();
            Some(Region {
                first: header.end_byte(),
                last: node.end_byte(),
                private,
                indices,
                ..Region::default()
            })
        })
        .collect()
}

fn independent_item(unit: &Unit, target: Node<'_>, region: &Region) -> bool {
    walk(target)
        .filter(|node| node.kind() == "call_expression")
        .any(|node| {
            let indices: BTreeSet<_> = arguments(unit, node)
                .into_iter()
                .map(|argument| unit.name(argument.node))
                .collect();
            !region.indices.is_empty() && region.indices.is_subset(&indices)
        })
}

fn static_storage(project: &Project<'_>, id: usize) -> bool {
    let symbol = &project.symbols[id];
    let scope = &project.scopes[symbol.scope];
    if matches!(scope.kind, "module" | "submodule" | "program")
        || symbol.attrs.contains("save")
        || matches!(
            declarator(symbol).kind(),
            "init_declarator" | "pointer_init_declarator"
        )
    {
        return true;
    }
    children(scope.node)
        .into_iter()
        .filter(|node| node.kind() == "save_statement")
        .any(|save| {
            let names = children(save);
            names.is_empty()
                || names
                    .into_iter()
                    .any(|node| project.units[scope.unit].name(node) == symbol.name)
        })
}

fn mutation_targets<'a>(procedure: &Procedure<'a, '_>, node: Node<'a>) -> Vec<Node<'a>> {
    match node.kind() {
        "assignment_statement" | "pointer_association_statement" => {
            if node
                .parent()
                .is_some_and(|parent| matches!(parent.kind(), "association" | "selector"))
            {
                return Vec::new();
            }
            node.child_by_field_name("left")
                .or_else(|| node.named_child(0))
                .into_iter()
                .collect()
        }
        "allocate_statement" => fields(node, "allocation")
            .into_iter()
            .map(|node| {
                if node.kind() == "sized_allocation" {
                    node.named_child(0).unwrap_or(node)
                } else {
                    node
                }
            })
            .collect(),
        "deallocate_statement" | "nullify_statement" => children(node)
            .into_iter()
            .filter(|part| {
                matches!(
                    part.kind(),
                    "identifier" | "derived_type_member_expression" | "call_expression"
                )
            })
            .collect(),
        "subroutine_call" => formal_arguments(procedure.project, procedure.scope, node)
            .into_iter()
            .filter_map(|(actual, id)| {
                let attrs = &procedure.project.symbols[id].attrs;
                (attrs.contains("intent(inout)") || attrs.contains("intent(out)")).then_some(actual)
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn side_effects(procedure: &Procedure<'_, '_>, nodes: &[Node<'_>]) -> Vec<Diagnostic> {
    let scope = &procedure.project.scopes[procedure.scope];
    let policy = &procedure.config.safety;
    let attributes = |name| {
        scope.attrs.contains(name)
            || scope
                .declaration
                .is_some_and(|id| procedure.project.scopes[id].attrs.contains(name))
    };
    let pure = policy.check_pure
        && (attributes("pure") || attributes("elemental"))
        && !attributes("impure");
    let (omp, protected) = if policy.check_openmp {
        openmp_regions(procedure.unit, nodes)
    } else {
        (Vec::new(), Vec::new())
    };
    let concurrent = if policy.check_do_concurrent {
        concurrent_regions(procedure.unit, nodes)
    } else {
        Vec::new()
    };
    if !pure && omp.is_empty() && concurrent.is_empty() {
        return Vec::new();
    }
    let mut diagnostics = Vec::new();
    for node in nodes.iter().copied() {
        for target in mutation_targets(procedure, node) {
            let Some(place) = procedure.place(target) else {
                continue;
            };
            let symbol = &procedure.project.symbols[place.root];
            let name = &symbol.name;
            let mut declaring = symbol.scope;
            while procedure.project.scopes[declaring].kind == "block" {
                let Some(parent) = procedure.project.scopes[declaring].parent else {
                    break;
                };
                declaring = parent;
            }
            // Separate module procedures bind the dummies and result declared in their interface.
            let inherited_signature = scope.declaration == Some(declaring)
                && scope.symbols.get(name) == Some(&place.root);
            if pure && declaring != procedure.scope && !inherited_signature {
                diagnostics.push(procedure.diagnostic(
                    "PRL006",
                    format!("PURE procedure modifies nonlocal storage '{name}'"),
                    target,
                    true,
                ));
                continue;
            }
            if !static_storage(procedure.project, place.root) {
                continue;
            }
            if concurrent.iter().any(|region| {
                region.contains(node)
                    && !region.private.contains(name)
                    && !independent_item(procedure.unit, target, region)
            }) {
                let certain = !place
                    .parts
                    .iter()
                    .any(|part| matches!(part, Part::Indices(_)));
                diagnostics.push(procedure.diagnostic(
                    "PRL006",
                    format!("DO CONCURRENT may modify shared static storage '{name}'"),
                    target,
                    certain,
                ));
                continue;
            }
            let open: Vec<_> = omp.iter().filter(|region| region.contains(node)).collect();
            if !open.is_empty()
                && !open
                    .iter()
                    .any(|region| region.private.contains(name) || region.private.contains("*"))
                && !protected.iter().any(|region| {
                    region.contains(node) && (region.protect || region.private.contains(name))
                })
            {
                diagnostics.push(procedure.diagnostic(
                    "PRL006",
                    format!("OpenMP parallel region may modify shared static storage '{name}'"),
                    target,
                    false,
                ));
            }
        }
    }
    diagnostics
}

pub fn check(project: &Project<'_>, config: &Config) -> Vec<Diagnostic> {
    let mut diagnostics = BTreeSet::new();
    for (id, scope) in project
        .scopes
        .iter()
        .enumerate()
        .filter(|(_, scope)| executable(scope))
    {
        let nodes = body_nodes(scope);
        let procedure = Procedure::new(project, config, id, &nodes);
        diagnostics.extend(side_effects(&procedure, &nodes));
        if !procedure.calls.is_empty() {
            let control = flow::build(scope, procedure.unit, |node| {
                config.safety.terminates(procedure.calls.get(&node.id()))
            });
            diagnostics.extend(procedure.diagnostics(&control));
        }
    }
    diagnostics.into_iter().collect()
}
