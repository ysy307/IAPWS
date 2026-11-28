use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;
use std::path::Path;

use serde_json::{Value, json};

use crate::config::{Config, matches};
use crate::diagnostic::{Diagnostic, Location};
use crate::semantic::{Project, Reference, Scope, Symbol, Use};
use crate::source::{child, children, walk};

pub struct DesignReport {
    pub diagnostics: Vec<Diagnostic>,
    pub metrics: Value,
    pub graph: String,
}

struct Method {
    name: String,
    components: BTreeSet<String>,
    calls: BTreeSet<String>,
    complexity: usize,
    location: Location,
}

struct TypeInfo {
    location: Location,
    base: Option<String>,
    abstract_type: bool,
    stateless: bool,
    methods: Vec<Method>,
    dependencies: BTreeSet<String>,
}

fn cyclomatic(scope: &Scope<'_>, config: &Config) -> usize {
    1 + children(scope.node)
        .into_iter()
        .filter(|node| !matches!(node.kind(), "internal_procedures" | "interface"))
        .flat_map(walk)
        .filter(|node| {
            config
                .metrics
                .complexity_nodes
                .iter()
                .any(|kind| kind == node.kind())
                && child(*node, "default").is_none()
        })
        .count()
}

fn scope_location(project: &Project<'_>, scope: usize) -> Location {
    let scope = &project.scopes[scope];
    project.units[scope.unit].location(scope.header)
}

fn ancestry(project: &Project<'_>, owner: usize) -> BTreeSet<usize> {
    let mut inherited = BTreeSet::new();
    let mut current = Some(owner);
    while let Some(id) = current {
        if !inherited.insert(id) {
            break;
        }
        current = project.resolve_parent(id);
    }
    inherited
}

fn passed_object(project: &Project<'_>, body: usize, binding: &Symbol<'_>) -> Option<usize> {
    if binding.attrs.contains("nopass") {
        return None;
    }
    let scope = &project.scopes[body];
    let declaration = &project.scopes[scope.declaration.unwrap_or(body)];
    let default = scope
        .parameters
        .first()
        .or_else(|| declaration.parameters.first())?;
    let name = binding
        .attrs
        .iter()
        .find_map(|attr| attr.strip_prefix("pass(")?.strip_suffix(')'));
    project.lookup(body, name.unwrap_or(default))
}

fn type_infos<'a>(
    project: &Project<'a>,
    config: &Config,
    references: &[Vec<&Reference<'a>>],
) -> (BTreeMap<String, TypeInfo>, BTreeSet<usize>) {
    let mut result = BTreeMap::new();
    let mut bound = BTreeSet::new();
    for &owner in &project.types {
        let ty = &project.symbols[owner];
        let state = ancestry(project, owner);
        let mut dependencies: BTreeSet<_> = ty
            .members
            .values()
            .filter_map(|&item| project.type_of(item))
            .collect();
        let mut bodies = Vec::new();
        let mut aliases: BTreeMap<usize, BTreeSet<String>> = BTreeMap::new();
        for &binding in ty.members.values() {
            let method = &project.symbols[binding];
            if method.kind != "binding" {
                continue;
            }
            let Some(target) = project.binding_target(binding) else {
                continue;
            };
            let implementation = project.implementation(target);
            let formal = project.symbols[target].procedure;
            for scope in [formal, implementation].into_iter().flatten() {
                for &item in project.scopes[scope].symbols.values() {
                    dependencies.extend(project.type_of(item));
                }
                for reference in &references[scope] {
                    dependencies.extend(project.type_of(reference.symbol));
                    if reference.call {
                        dependencies.extend(project.symbols[reference.symbol].owner);
                    }
                }
            }
            if let Some(body) = implementation.filter(|_| !method.attrs.contains("deferred")) {
                bound.insert(body);
                bodies.push((binding, body));
                aliases
                    .entry(target)
                    .or_default()
                    .insert(method.qualified.clone());
                let scope = &project.scopes[body];
                if let Some(other) = scope
                    .parent
                    .and_then(|p| project.scopes[p].symbols.get(&scope.name))
                {
                    aliases
                        .entry(*other)
                        .or_default()
                        .insert(method.qualified.clone());
                }
            }
        }
        let mut methods = Vec::new();
        for (binding, body) in bodies {
            let method = &project.symbols[binding];
            let object = passed_object(project, body, method);
            let mut components = BTreeSet::new();
            let mut calls = BTreeSet::new();
            for reference in &references[body] {
                let target = &project.symbols[reference.symbol];
                let object_member = object.is_some()
                    && reference.root_receiver == object
                    && target.owner.is_some_and(|parent| state.contains(&parent));
                let accessed_through_parent = reference
                    .receiver
                    .is_some_and(|receiver| project.symbols[receiver].kind == "parent_component");
                if object_member
                    && matches!(target.kind, "variable" | "procedure")
                    && (reference.receiver == object || accessed_through_parent)
                {
                    components.insert(target.qualified.clone());
                }
                if reference.call {
                    if let Some(names) = aliases.get(&reference.symbol) {
                        calls.extend(names.iter().cloned());
                        continue;
                    }
                    let actual = if object_member && reference.receiver == object {
                        project
                            .member(owner, &target.name)
                            .map_or(target, |id| &project.symbols[id])
                    } else {
                        target
                    };
                    calls.insert(actual.qualified.clone());
                }
            }
            methods.push(Method {
                name: method.qualified.clone(),
                components,
                calls,
                complexity: cyclomatic(&project.scopes[body], config),
                location: scope_location(project, body),
            });
        }
        methods.sort_unstable_by(|a, b| a.name.cmp(&b.name));
        let base = project
            .resolve_parent(owner)
            .map(|id| project.symbols[id].qualified.clone())
            .or_else(|| ty.base_name.clone());
        let stateless = !state.iter().any(|&id| {
            project.symbols[id].members.values().any(|&member| {
                let member = &project.symbols[member];
                matches!(member.kind, "variable" | "procedure")
                    && !member.attrs.contains("kind")
                    && !member.attrs.contains("len")
            })
        });
        let mut dependencies: BTreeSet<_> = dependencies
            .into_iter()
            .filter(|&id| id != owner)
            .map(|id| project.symbols[id].qualified.clone())
            .collect();
        dependencies.extend(base.iter().filter(|name| *name != &ty.qualified).cloned());
        result
            .entry(ty.qualified.clone())
            .or_insert_with(|| TypeInfo {
                location: project.units[project.scopes[ty.scope].unit].location(ty.node),
                base,
                abstract_type: ty.attrs.contains("abstract"),
                stateless,
                methods,
                dependencies,
            });
    }
    (result, bound)
}

fn inheritance<'a>(
    types: &'a BTreeMap<String, TypeInfo>,
    diagnostics: &mut Vec<Diagnostic>,
) -> (BTreeMap<&'a str, Option<usize>>, BTreeMap<&'a str, usize>) {
    let mut depths = BTreeMap::new();
    let mut children = BTreeMap::new();
    for info in types.values() {
        if let Some(base) = info.base.as_deref() {
            *children.entry(base).or_insert(0) += 1;
        }
    }
    for start in types.keys() {
        let mut chain = Vec::new();
        let mut positions = BTreeMap::new();
        let mut next = start.as_str();
        let mut depth = loop {
            if let Some(&depth) = depths.get(next) {
                break depth;
            }
            let info = &types[next];
            if let Some(begin) = positions.insert(next, chain.len()) {
                let mut cycle: Vec<&str> = chain[begin..].to_vec();
                cycle.push(next);
                diagnostics.push(
                    Diagnostic::new(
                        "DIT001",
                        format!("Circular type extension involving {}", cycle.join(" -> ")),
                        info.location.clone(),
                    )
                    .with_subject(next),
                );
                break None;
            }
            let Some(parent) = info.base.as_deref() else {
                depths.insert(next, Some(0));
                break Some(0);
            };
            chain.push(next);
            if !types.contains_key(parent) {
                diagnostics.push(
                    Diagnostic::new(
                        "DIT002",
                        format!("Cannot determine the inheritance depth of base type {parent}"),
                        info.location.clone(),
                    )
                    .with_subject(next),
                );
                break None;
            }
            next = parent;
        };
        for item in chain.into_iter().rev() {
            depth = depth.map(|value| value + 1);
            depths.insert(item, depth);
        }
    }
    (depths, children)
}

fn set_root(parents: &mut [usize], mut id: usize) -> usize {
    while parents[id] != id {
        parents[id] = parents[parents[id]];
        id = parents[id];
    }
    id
}

fn method_components(methods: &[Method]) -> Vec<Vec<&str>> {
    let mut parents: Vec<_> = (0..methods.len()).collect();
    let indices: BTreeMap<_, _> = methods
        .iter()
        .enumerate()
        .map(|(i, m)| (m.name.as_str(), i))
        .collect();
    let mut first_access = BTreeMap::new();
    for (id, method) in methods.iter().enumerate() {
        for component in &method.components {
            let other = *first_access.entry(component.as_str()).or_insert(id);
            let found = set_root(&mut parents, other);
            parents[found] = set_root(&mut parents, id);
        }
        for call in &method.calls {
            if let Some(&other) = indices.get(call.as_str()) {
                let found = set_root(&mut parents, other);
                parents[found] = set_root(&mut parents, id);
            }
        }
    }
    let mut groups: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for (id, method) in methods.iter().enumerate() {
        groups
            .entry(set_root(&mut parents, id))
            .or_default()
            .push(method.name.as_str());
    }
    let mut groups: Vec<_> = groups.into_values().collect();
    groups.sort_unstable();
    groups
}

fn limit(
    diagnostics: &mut Vec<Diagnostic>,
    code: &str,
    label: &str,
    name: &str,
    value: usize,
    max: usize,
    location: &Location,
) {
    if value > max {
        let diagnostic = Diagnostic::new(
            code,
            format!("{name} has {label} {value} (limit {max})"),
            location.clone(),
        )
        .with_subject(name);
        diagnostics.push(if matches!(code, "CBO001" | "WMC001" | "CC001") {
            diagnostic.warning()
        } else {
            diagnostic
        });
    }
}

fn type_metrics(
    project: &Project<'_>,
    config: &Config,
    diagnostics: &mut Vec<Diagnostic>,
) -> BTreeMap<String, Value> {
    let mut references = vec![Vec::new(); project.scopes.len()];
    for item in &project.references {
        let mut scope = item.scope;
        while let ("block", Some(parent)) =
            (project.scopes[scope].kind, project.scopes[scope].parent)
        {
            scope = parent;
        }
        references[scope].push(item);
    }
    let (types, bound) = type_infos(project, config, &references);
    let (depths, children) = inheritance(&types, diagnostics);
    let mut metrics = BTreeMap::new();
    let max = &config.limits;
    for (name, info) in &types {
        let subject_limits = max.subjects.get(name);
        let depth = depths[name.as_str()];
        if let Some(depth) = depth {
            limit(
                diagnostics,
                "DIT003",
                "inheritance depth",
                name,
                depth,
                max.dit,
                &info.location,
            );
        }
        let components = method_components(&info.methods);
        let lcom4 = components.len();
        if !(info.abstract_type && config.metrics.exclude_abstract
            || info.stateless && config.metrics.exclude_stateless)
        {
            limit(
                diagnostics,
                "COH001",
                "LCOM4",
                name,
                lcom4,
                subject_limits
                    .and_then(|limits| limits.lcom4)
                    .unwrap_or(max.lcom4),
                &info.location,
            );
        }
        let cbo = info.dependencies.len();
        let wmc = info.methods.iter().map(|method| method.complexity).sum();
        let mut largest = 0;
        let mut responses = BTreeSet::new();
        let mut calls = BTreeMap::new();
        let mut accesses = BTreeMap::new();
        for method in &info.methods {
            largest = largest.max(method.complexity);
            responses.insert(&method.name);
            responses.extend(&method.calls);
            accesses.insert(&method.name, &method.components);
            calls.insert(&method.name, &method.calls);
            limit(
                diagnostics,
                "CC001",
                "cyclomatic complexity",
                &method.name,
                method.complexity,
                max.complexity_for(&method.name),
                &method.location,
            );
        }
        limit(
            diagnostics,
            "CBO001",
            "CBO",
            name,
            cbo,
            subject_limits.and_then(|limits| limits.cbo).unwrap_or(max.cbo),
            &info.location,
        );
        limit(
            diagnostics,
            "WMC001",
            "WMC",
            name,
            wmc,
            subject_limits.and_then(|limits| limits.wmc).unwrap_or(max.wmc),
            &info.location,
        );
        metrics.insert(name.clone(), json!({
            "inheritance": {"base": info.base, "dit": depth, "noc": children.get(name.as_str()).copied().unwrap_or(0)},
            "cohesion": {"lcom4": lcom4, "components": components, "calls": calls, "component_accesses": accesses},
            "coupling": {"cbo": cbo, "wmc": wmc, "rfc": responses.len(), "max_complexity": largest, "dependencies": info.dependencies},
        }));
    }
    for (ind, scope) in project.scopes.iter().enumerate() {
        if matches!(scope.kind, "subroutine" | "function" | "module_procedure")
            && !scope.interface_body
            && !bound.contains(&ind)
        {
            limit(
                diagnostics,
                "CC001",
                "cyclomatic complexity",
                &scope.qualified,
                cyclomatic(scope, config),
                max.complexity_for(&scope.qualified),
                &scope_location(project, ind),
            );
        }
    }
    metrics
}

#[derive(Default)]
struct ModuleInfo {
    uses: Vec<Use>,
    external: Vec<Use>,
}

fn public_units(project: &Project<'_>, config: &Config) -> BTreeMap<usize, ModuleInfo> {
    let mut units = BTreeMap::new();
    let mut roots = Vec::with_capacity(project.scopes.len());
    for (id, scope) in project.scopes.iter().enumerate() {
        let mut root = id;
        while let Some(parent) = project.scopes[root].parent {
            root = parent;
        }
        roots.push(root);
        let module = units.entry(root).or_insert_with(ModuleInfo::default);
        module.uses.extend(scope.uses.iter().cloned());
        if scope.kind == "submodule"
            && let Some(ancestor) = &scope.ancestor
        {
            let name = scope
                .parent_name
                .as_ref()
                .map_or_else(|| ancestor.clone(), |parent| format!("{ancestor}:{parent}"));
            module
                .uses
                .push(dependency(name, scope_location(project, id), false));
        }
    }
    let external =
        project
            .external_imports
            .iter()
            .cloned()
            .chain(project.references.iter().filter_map(|reference| {
                project.symbols[reference.symbol].external.as_ref()?;
                let scope = &project.scopes[reference.scope];
                Some((
                    reference.scope,
                    reference.symbol,
                    project.units[scope.unit].location(reference.node),
                ))
            }));
    for (scope, id, site) in external {
        let symbol = &project.symbols[id];
        let Some(origin) = symbol.external.as_deref() else {
            continue;
        };
        let name = if origin.starts_with("petsc")
            && symbol.name.starts_with("mpi_")
            && config.externals.contains_key("mpi")
        {
            "mpi"
        } else {
            origin
        };
        if let Some(module) = config.externals.get(name) {
            units
                .get_mut(&roots[scope])
                .expect("each scope has a module root")
                .external
                .push(dependency(name.to_owned(), site, module.intrinsic));
        }
    }
    for info in units.values_mut() {
        info.external
            .sort_unstable_by(|a, b| compare_locations(&a.location, &b.location));
    }
    units
}

fn dependency(module: String, location: Location, intrinsic: bool) -> Use {
    Use {
        module,
        location,
        intrinsic,
        non_intrinsic: false,
        only: false,
        imports: BTreeMap::new(),
        imported_at: BTreeMap::new(),
    }
}

fn compare_locations(a: &Location, b: &Location) -> std::cmp::Ordering {
    (&a.path, a.line, a.column).cmp(&(&b.path, b.line, b.column))
}

fn unit_name(scope: &Scope<'_>) -> String {
    if matches!(scope.kind, "module" | "submodule") {
        return scope.name.to_ascii_lowercase();
    }
    format!("{}:{}", scope.kind, scope.name.to_ascii_lowercase())
}

fn path_in_root(path: &Path, root: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
        .to_lowercase()
}

fn cyclic_components(graph: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut reverse = vec![Vec::new(); graph.len()];
    for (from, outgoing) in graph.iter().enumerate() {
        for &to in outgoing {
            reverse[to].push(from);
        }
    }
    let mut visited = vec![false; graph.len()];
    let mut finished = Vec::new();
    for start in 0..graph.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![(start, 0)];
        while let Some((node, index)) = stack.last_mut() {
            if let Some(&next) = graph[*node].get(*index) {
                *index += 1;
                if !visited[next] {
                    visited[next] = true;
                    stack.push((next, 0));
                }
            } else {
                finished.push(*node);
                stack.pop();
            }
        }
    }
    visited.fill(false);
    let mut components = Vec::new();
    for start in finished.into_iter().rev() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![start];
        let mut component = Vec::new();
        while let Some(node) = stack.pop() {
            component.push(node);
            for &next in &reverse[node] {
                if !visited[next] {
                    visited[next] = true;
                    stack.push(next);
                }
            }
        }
        if component.len() > 1 || graph[start].contains(&start) {
            component.sort_unstable();
            components.push(component);
        }
    }
    components.sort_unstable();
    components
}

type Graph = BTreeMap<String, BTreeSet<String>>;

fn dependency_metrics(
    project: &Project<'_>,
    config: &Config,
    root: &Path,
    diagnostics: &mut Vec<Diagnostic>,
) -> (BTreeMap<String, Value>, Graph, String) {
    let units = public_units(project, config);
    let mut sorted: Vec<_> = units
        .keys()
        .map(|&id| {
            (
                unit_name(&project.scopes[id]),
                scope_location(project, id),
                id,
            )
        })
        .collect();
    sorted.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| compare_locations(&a.1, &b.1)));
    let mut definitions: BTreeMap<String, (Location, usize)> = BTreeMap::new();
    for (name, location, id) in sorted {
        if let Some((first, _)) = definitions.get(&name) {
            diagnostics.push(
                Diagnostic::new(
                    "ARCH002",
                    format!(
                        "Duplicate module {name} (first definition at {}:{})",
                        first.path.display(),
                        first.line
                    ),
                    location,
                )
                .with_subject(&name),
            );
        } else {
            definitions.insert(name, (location, id));
        }
    }

    let paths: BTreeMap<_, _> = definitions
        .iter()
        .map(|(name, (site, _))| (name.as_str(), path_in_root(&site.path, root)))
        .collect();
    let mut layers = BTreeMap::new();
    for (name, (location, _)) in &definitions {
        let path = &paths[name.as_str()];
        let matching: Vec<_> = config
            .layers
            .iter()
            .filter(|(_, layer)| layer.paths.iter().any(|pattern| matches(pattern, path)))
            .map(|(item, _)| item.as_str())
            .collect();
        if let [single] = matching.as_slice() {
            layers.insert(name.as_str(), Some(*single));
        } else {
            layers.insert(name.as_str(), None);
            if !config.layers.is_empty() {
                let message = if matching.is_empty() {
                    format!("No architecture layer matches {path}")
                } else {
                    format!(
                        "Multiple architecture layers ({}) match {path}",
                        matching.join(", ")
                    )
                };
                diagnostics
                    .push(Diagnostic::new("ARCH001", message, location.clone()).with_subject(name));
            }
        }
    }

    let mut graph = Graph::new();
    let mut sites = BTreeMap::new();
    for (origin, (_, id)) in &definitions {
        let layer = layers[origin.as_str()];
        let mut targets = BTreeSet::new();
        let mut external_uses = BTreeSet::new();
        let mut internal_uses = BTreeSet::new();
        let info = &units[id];
        for (usage, reexport) in info
            .uses
            .iter()
            .map(|item| (item, false))
            .chain(info.external.iter().map(|item| (item, true)))
        {
            let target = usage.module.to_ascii_lowercase();
            let subject = format!("{origin}::{target}");
            if definitions.contains_key(&target) && !usage.intrinsic && !reexport {
                if !internal_uses.insert(target.clone()) {
                    continue;
                }
                targets.insert(target.clone());
                sites.insert((origin.clone(), target.clone()), &usage.location);
                if let (Some(from), Some(to)) = (layer, layers[target.as_str()])
                    && !config.layers[from]
                        .may_use
                        .iter()
                        .any(|allowed| allowed == "*" || allowed == to)
                {
                    diagnostics.push(
                        Diagnostic::new(
                            "ARCH004",
                            format!("Layer {from} may not use {target} in layer {to}"),
                            usage.location.clone(),
                        )
                        .with_subject(&subject),
                    );
                }
                if config.internals.iter().any(|item| {
                    matches(&item.module, &target)
                        && !item
                            .layers
                            .iter()
                            .any(|allowed| allowed == "*" || Some(allowed.as_str()) == layer)
                }) {
                    let from = layer.unwrap_or("<unassigned>");
                    diagnostics.push(
                        Diagnostic::new(
                            "ARCH005",
                            format!("Layer {from} may not access internal module {target}"),
                            usage.location.clone(),
                        )
                        .with_subject(&subject),
                    );
                }
                continue;
            }
            let external = config
                .externals
                .get(&target)
                .filter(|ext| !usage.intrinsic || ext.intrinsic)
                .filter(|ext| !usage.non_intrinsic || !ext.intrinsic);
            let Some(external) = external else {
                diagnostics.push(
                    Diagnostic::new(
                        "ARCH003",
                        format!("Cannot resolve module {target}"),
                        usage.location.clone(),
                    )
                    .with_subject(&subject),
                );
                continue;
            };
            if !external_uses.insert(target.clone()) {
                continue;
            }
            targets.insert(format!("external::{target}"));
            let path = &paths[origin.as_str()];
            let allowed_path = external.paths.is_empty()
                || external.paths.iter().any(|pattern| matches(pattern, path));
            let allowed_layer = external.layers.is_empty()
                || external
                    .layers
                    .iter()
                    .any(|allowed| allowed == "*" || Some(allowed.as_str()) == layer);
            if !allowed_path || !allowed_layer {
                diagnostics.push(
                    Diagnostic::new(
                        "ARCH006",
                        format!("{path} may not access external module {target}"),
                        usage.location.clone(),
                    )
                    .with_subject(&subject),
                );
            }
        }
        graph.insert(origin.clone(), targets);
    }

    let names: Vec<_> = graph.keys().collect();
    let indices: BTreeMap<_, _> = names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let internal: Vec<Vec<_>> = graph
        .values()
        .map(|outgoing| {
            outgoing
                .iter()
                .filter_map(|item| indices.get(item.as_str()).copied())
                .collect()
        })
        .collect();
    for group in cyclic_components(&internal) {
        let first = group[0];
        let members: Vec<_> = group.iter().map(|&id| names[id].as_str()).collect();
        let neighbor = internal[first]
            .iter()
            .find(|next| group.binary_search(next).is_ok())
            .expect("cyclic component has outgoing edge");
        diagnostics.push(
            Diagnostic::new(
                "ARCH007",
                format!("Cyclic module dependencies among {}", members.join(", ")),
                sites[&(names[first].clone(), names[*neighbor].clone())].clone(),
            )
            .with_subject(names[first]),
        );
    }

    let mut inbound = vec![0; names.len()];
    for links in &internal {
        for &target in links {
            inbound[target] += 1;
        }
    }
    let mut metrics = BTreeMap::new();
    for (ind, &name) in names.iter().enumerate() {
        let fan_in = inbound[ind];
        let fan_out = if config.metrics.external_fan_out {
            graph[name].len()
        } else {
            internal[ind].len()
        };
        let instability = if fan_in + fan_out == 0 {
            0.0
        } else {
            fan_out as f64 / (fan_in + fan_out) as f64
        };
        metrics.insert(name.clone(), json!({"layer": layers[name.as_str()], "fan_in": fan_in, "fan_out": fan_out, "instability": instability}));
    }
    let dot = graphviz(&graph, &layers);
    (metrics, graph, dot)
}

fn quote(label: &str) -> String {
    serde_json::to_string(label).expect("a string can always be quoted")
}

fn graphviz(graph: &Graph, layers: &BTreeMap<&str, Option<&str>>) -> String {
    let mut text = String::from("digraph fortran_modules {\n    rankdir=LR;\n");
    let mut grouped: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for (&name, &layer) in layers {
        if let Some(layer) = layer {
            grouped.entry(layer).or_default().push(name);
        } else {
            let _ = writeln!(text, "    {};", quote(name));
        }
    }
    for (layer, names) in grouped {
        let _ = writeln!(
            text,
            "    subgraph {} {{",
            quote(&format!("cluster_{layer}"))
        );
        let _ = writeln!(text, "        label={};", quote(layer));
        for name in names {
            let _ = writeln!(text, "        {};", quote(name));
        }
        text.push_str("    }\n");
    }
    let externals: BTreeSet<_> = graph
        .values()
        .flatten()
        .filter(|item| !graph.contains_key(*item))
        .collect();
    for name in externals {
        let label = name.strip_prefix("external::").unwrap_or(name);
        let _ = writeln!(
            text,
            "    {} [label={},shape=box,style=dashed];",
            quote(name),
            quote(label)
        );
    }
    for (from, outgoing) in graph {
        for to in outgoing {
            let _ = writeln!(text, "    {} -> {};", quote(from), quote(to));
        }
    }
    text.push_str("}\n");
    text
}

pub fn analyze(project: &Project<'_>, config: &Config, root: &Path) -> DesignReport {
    let mut diagnostics = Vec::new();
    let types = type_metrics(project, config, &mut diagnostics);
    let (modules, edges, graph) = dependency_metrics(project, config, root, &mut diagnostics);
    DesignReport {
        diagnostics,
        metrics: json!({"types": types, "modules": modules, "module_dependencies": edges}),
        graph,
    }
}
