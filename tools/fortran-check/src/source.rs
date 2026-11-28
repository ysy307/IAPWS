use std::path::PathBuf;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use tree_sitter::{Node, Parser, Tree};

use crate::diagnostic::{Diagnostic, Location};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    pub path: PathBuf,
    pub text: String,
    pub line_map: Vec<Option<(PathBuf, usize)>>,
}

impl Source {
    pub fn new(path: impl Into<PathBuf>, text: impl Into<String>) -> Self {
        Self {
            path: path.into(),
            text: text.into(),
            line_map: Vec::new(),
        }
    }

    pub fn location(&self, row: usize, column: usize) -> Location {
        match row
            .checked_sub(1)
            .and_then(|index| self.line_map.get(index))
            .and_then(Option::as_ref)
        {
            Some((path, line)) => Location::new(path, *line, column),
            None => Location::new(&self.path, row, column),
        }
    }
}

#[derive(Debug)]
pub struct Unit {
    pub source: Source,
    pub tree: Tree,
}

impl Unit {
    pub fn parse(source: Source) -> Result<Self> {
        let mut parser = Parser::new();
        parser
            .set_language(&tree_sitter_fortran::LANGUAGE.into())
            .context("could not load the Fortran grammar")?;
        let tree = parser
            .parse(&source.text, None)
            .with_context(|| format!("no parse tree for {}", source.path.display()))?;
        Ok(Self { source, tree })
    }

    pub fn text(&self, node: Node<'_>) -> &str {
        self.source.text.get(node.byte_range()).unwrap_or_default()
    }

    pub fn name(&self, node: Node<'_>) -> String {
        self.text(node)
            .chars()
            .filter(|part| !part.is_whitespace())
            .map(|part| part.to_ascii_lowercase())
            .collect()
    }

    pub fn location(&self, node: Node<'_>) -> Location {
        let pos = node.start_position();
        self.source.location(pos.row + 1, pos.column + 1)
    }

    pub fn diagnostics(&self) -> Vec<Diagnostic> {
        let mut pending = vec![self.tree.root_node()];
        let mut diagnostics = Vec::new();
        while let Some(node) = pending.pop() {
            if node.is_missing() {
                diagnostics.push(Diagnostic::new(
                    "PAR001",
                    format!("Incomplete Fortran syntax (expected {})", node.kind()),
                    self.location(node),
                ));
            } else if node.is_error() {
                diagnostics.push(Diagnostic::new(
                    "PAR001",
                    "Cannot parse Fortran syntax",
                    self.location(node),
                ));
            } else if node.kind().starts_with("preproc_") || node.kind() == "include_statement" {
                diagnostics.push(Diagnostic::new(
                    "PAR002",
                    "Unexpanded include or preprocessor directive",
                    self.location(node),
                ));
            } else if node.has_error() {
                let mut cursor = node.walk();
                let children: Vec<_> = node.children(&mut cursor).collect();
                pending.extend(children.into_iter().rev());
            } else if node.is_named() {
                pending.extend(children(node).into_iter().rev());
            }
        }
        diagnostics
    }
}

pub fn parse(source: Source) -> Result<Unit> {
    Unit::parse(source)
}

pub fn children(node: Node<'_>) -> Vec<Node<'_>> {
    let mut cursor = node.walk();
    node.named_children(&mut cursor)
        .filter(|child| child.kind() != "comment")
        .collect()
}

pub fn child<'a>(node: Node<'a>, kind: &str) -> Option<Node<'a>> {
    let mut cursor = node.walk();
    node.named_children(&mut cursor)
        .find(|item| item.kind() == kind)
}

pub fn fields<'a>(node: Node<'a>, field: &str) -> Vec<Node<'a>> {
    let mut cursor = node.walk();
    node.children_by_field_name(field, &mut cursor).collect()
}

pub fn walk(node: Node<'_>) -> impl Iterator<Item = Node<'_>> {
    let mut stack = vec![node];
    std::iter::from_fn(move || {
        let current = stack.pop()?;
        let mut cursor = current.walk();
        let children: Vec<_> = current.named_children(&mut cursor).collect();
        stack.extend(children.into_iter().rev());
        Some(current)
    })
}
