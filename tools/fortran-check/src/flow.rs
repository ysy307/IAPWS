use std::collections::{BTreeSet, VecDeque};

use tree_sitter::Node;

use crate::{
    semantic::Scope,
    source::{Unit, child, children},
};

#[derive(Debug, Clone, Copy)]
pub struct Edge<'a> {
    pub target: usize,
    pub condition: Option<Node<'a>>,
    pub truth: bool,
}

impl Edge<'_> {
    fn direct(target: usize) -> Self {
        Self {
            target,
            condition: None,
            truth: true,
        }
    }
}

#[derive(Debug)]
pub struct Step<'a> {
    pub node: Option<Node<'a>>,
    pub edges: Vec<Edge<'a>>,
}

#[derive(Debug)]
pub struct Flow<'a> {
    pub entry: usize,
    pub exit: usize,
    pub steps: Vec<Step<'a>>,
}

struct Construct {
    name: Option<String>,
    after: usize,
    repeat: Option<usize>,
}

pub fn unwrap(mut node: Node<'_>) -> Node<'_> {
    while node.kind() == "parenthesized_expression" {
        let Some(inner) = children(node)
            .into_iter()
            .find(|part| part.kind() != "comment")
        else {
            break;
        };
        node = inner;
    }
    node
}

struct Builder<'a, 'b, F> {
    unit: &'b Unit,
    steps: Vec<Step<'a>>,
    nesting: Vec<Construct>,
    terminates: F,
}

impl<'a, F: Fn(Node<'a>) -> bool> Builder<'a, '_, F> {
    fn step(&mut self, node: Option<Node<'a>>, edges: Vec<Edge<'a>>) -> usize {
        let index = self.steps.len();
        self.steps.push(Step { node, edges });
        index
    }

    fn sequence(&mut self, nodes: Vec<Node<'a>>, after: usize) -> usize {
        nodes.into_iter().rev().fold(after, |next, node| {
            let kind = node.kind();
            if kind.starts_with("end_")
                || matches!(
                    kind,
                    "comment"
                        | "contains_statement"
                        | "derived_type_definition"
                        | "implicit_statement"
                        | "import_statement"
                        | "interface"
                        | "internal_procedures"
                        | "private_statement"
                        | "public_statement"
                        | "statement_label"
                        | "use_statement"
                        | "block_label_start_expression"
                        | "association_list"
                )
            {
                return next;
            }
            self.statement(node, next)
        })
    }

    fn edges(&self, condition: Node<'a>, when_true: usize, when_false: usize) -> Vec<Edge<'a>> {
        let branch = |truth| Edge {
            target: if truth { when_true } else { when_false },
            condition: Some(condition),
            truth,
        };
        let value = unwrap(condition);
        if value.kind() == "boolean_literal" {
            return vec![branch(self.unit.name(value).starts_with(".true."))];
        }
        vec![branch(true), branch(false)]
    }

    fn branch(&mut self, condition: Node<'a>, when_true: usize, when_false: usize) -> usize {
        self.step(
            Some(condition),
            self.edges(condition, when_true, when_false),
        )
    }

    fn enter(&mut self, node: Node<'a>, after: usize, repeat: Option<usize>) -> bool {
        let name = child(node, "block_label_start_expression")
            .map(|label| self.unit.name(label).trim_end_matches(':').to_owned());
        let pushed = name.is_some() || repeat.is_some();
        if pushed {
            self.nesting.push(Construct {
                name,
                after,
                repeat,
            });
        }
        pushed
    }

    fn leave(&mut self, pushed: bool) {
        if pushed {
            self.nesting.pop();
        }
    }

    fn conditional(&mut self, node: Node<'a>, after: usize) -> usize {
        let Some(condition) = child(node, "parenthesized_expression") else {
            return self.step(Some(node), vec![Edge::direct(after)]);
        };
        let is_clause = |part: &Node| {
            matches!(
                part.kind(),
                "elseif_clause" | "else_clause" | "elsewhere_clause"
            )
        };
        let pushed = self.enter(node, after, None);
        let parts = children(node);
        let mut fallback = after;
        for clause in parts.iter().rev().filter(|part| is_clause(part)) {
            let condition = child(*clause, "parenthesized_expression");
            let body = children(*clause)
                .into_iter()
                .filter(|part| Some(*part) != condition)
                .collect();
            let start = self.sequence(body, after);
            fallback = condition.map_or(start, |test| self.branch(test, start, fallback));
        }
        let body = parts
            .into_iter()
            .filter(|part| *part != condition && !is_clause(part))
            .collect();
        let when_true = self.sequence(body, after);
        self.leave(pushed);
        self.branch(condition, when_true, fallback)
    }

    fn loop_statement(&mut self, node: Node<'a>, after: usize) -> usize {
        let Some(header) = child(node, "do_statement") else {
            return self.step(Some(node), vec![Edge::direct(after)]);
        };
        let test = self.step(Some(header), Vec::new());
        let pushed = self.enter(node, after, Some(test));
        let body = children(node)
            .into_iter()
            .filter(|part| *part != header)
            .collect();
        let back = self.sequence(body, test);
        self.leave(pushed);
        let condition = child(header, "while_statement")
            .and_then(|while_part| child(while_part, "parenthesized_expression"));
        let mut edges = vec![Edge::direct(back)];
        if let Some(condition) = condition {
            edges = self.edges(condition, back, after);
        } else if children(header).iter().any(|part| {
            matches!(
                part.kind(),
                "loop_control_expression" | "concurrent_statement"
            )
        }) {
            edges.push(Edge::direct(after));
        }
        self.steps[test] = Step {
            node: Some(condition.unwrap_or(header)),
            edges,
        };
        test
    }

    fn select(&mut self, node: Node<'a>, after: usize, kind: &str) -> usize {
        let pushed = self.enter(node, after, None);
        let cases: Vec<_> = children(node)
            .into_iter()
            .filter(|part| part.kind() == kind)
            .collect();
        let mut targets = BTreeSet::new();
        for case in &cases {
            let type_selector = case.child_by_field_name("type");
            let body = children(*case)
                .into_iter()
                .filter(|part| {
                    Some(*part) != type_selector
                        && !matches!(part.kind(), "case_value_range_list" | "default")
                })
                .collect();
            targets.insert(self.sequence(body, after));
        }
        if !cases.iter().any(|case| child(*case, "default").is_some()) {
            targets.insert(after);
        }
        self.leave(pushed);
        self.step(
            child(node, "selector"),
            targets.into_iter().map(Edge::direct).collect(),
        )
    }

    fn statement(&mut self, node: Node<'a>, after: usize) -> usize {
        if (self.terminates)(node) {
            return self.step(Some(node), Vec::new());
        }
        match node.kind() {
            "if_statement" | "where_statement" => self.conditional(node, after),
            "do_loop" => self.loop_statement(node, after),
            "select_case_statement" => self.select(node, after, "case_statement"),
            "select_rank_statement" => self.select(node, after, "rank_statement"),
            "select_type_statement" => self.select(node, after, "type_statement"),
            "block_construct" | "associate_statement" => {
                let pushed = self.enter(node, after, None);
                let start = self.sequence(children(node), after);
                self.leave(pushed);
                start
            }
            "stop_statement" => self.step(Some(node), Vec::new()),
            "keyword_statement" => {
                let keyword = node.child(0).map(|part| part.kind());
                let name = child(node, "identifier").map(|part| self.unit.name(part));
                let matching = self.nesting.iter().rev().find(|construct| {
                    name.as_ref().map_or(construct.repeat.is_some(), |label| {
                        construct.name.as_ref() == Some(label)
                    })
                });
                let target = match keyword {
                    Some("return") => 0,
                    Some("exit") => matching.map_or(after, |construct| construct.after),
                    Some("cycle") => matching
                        .and_then(|construct| construct.repeat)
                        .unwrap_or(after),
                    _ => after,
                };
                self.step(Some(node), vec![Edge::direct(target)])
            }
            _ => self.step(Some(node), vec![Edge::direct(after)]),
        }
    }
}

pub fn build<'a>(
    scope: &Scope<'a>,
    unit: &Unit,
    terminates: impl Fn(Node<'a>) -> bool,
) -> Flow<'a> {
    let mut builder = Builder {
        unit,
        steps: Vec::new(),
        nesting: Vec::new(),
        terminates,
    };
    let exit = builder.step(None, Vec::new());
    let body = children(scope.node)
        .into_iter()
        .filter(|part| *part != scope.header)
        .collect();
    let start = builder.sequence(body, exit);
    let entry = builder.step(None, vec![Edge::direct(start)]);
    Flow {
        entry,
        exit,
        steps: builder.steps,
    }
}

pub fn solve<'a, State: Clone + PartialEq>(
    flow: &Flow<'a>,
    initial: State,
    transfer: impl Fn(&Step<'a>, &State) -> State,
    join: impl Fn(&State, &State) -> State,
    refine: impl Fn(&Edge<'a>, &State) -> Option<State>,
) -> Vec<Option<State>> {
    let mut facts = vec![None; flow.steps.len()];
    let mut queued = vec![false; flow.steps.len()];
    let mut pending = VecDeque::from([flow.entry]);
    facts[flow.entry] = Some(initial);
    queued[flow.entry] = true;
    while let Some(point) = pending.pop_front() {
        queued[point] = false;
        let step = &flow.steps[point];
        let after = transfer(step, facts[point].as_ref().expect("queued flow has input"));
        for edge in &step.edges {
            let Some(value) = refine(edge, &after) else {
                continue;
            };
            let previous = &mut facts[edge.target];
            let joined = previous
                .as_ref()
                .map_or_else(|| value.clone(), |old| join(old, &value));
            if previous.as_ref() == Some(&joined) {
                continue;
            }
            *previous = Some(joined);
            if !queued[edge.target] {
                queued[edge.target] = true;
                pending.push_back(edge.target);
            }
        }
    }
    facts
}
