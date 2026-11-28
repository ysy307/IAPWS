use std::{
    fmt,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Off,
    Note,
    Warning,
    #[default]
    Error,
}

impl fmt::Display for Severity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Off => "off",
            Self::Note => "note",
            Self::Warning => "warning",
            Self::Error => "error",
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Location {
    pub path: PathBuf,
    pub line: usize,
    pub column: usize,
}

impl Location {
    pub fn new(path: impl Into<PathBuf>, line: usize, column: usize) -> Self {
        Self {
            path: path.into(),
            line: line.max(1),
            column: column.max(1),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Diagnostic {
    pub location: Location,
    pub code: String,
    pub severity: Severity,
    pub message: String,
    pub subject: String,
}

impl Diagnostic {
    pub fn new(code: impl Into<String>, message: impl Into<String>, location: Location) -> Self {
        Self {
            location,
            code: code.into(),
            severity: Severity::Error,
            message: message.into(),
            subject: String::new(),
        }
    }

    pub fn with_subject(mut self, subject: impl Into<String>) -> Self {
        self.subject = subject.into();
        self
    }

    pub fn warning(mut self) -> Self {
        self.severity = Severity::Warning;
        self
    }

    pub fn gcc(&self) -> String {
        format!(
            "{}:{}:{}: {}: {} [{}]",
            self.location.path.display(),
            self.location.line,
            self.location.column,
            self.severity,
            self.message.replace(['\r', '\n'], " "),
            self.code
        )
    }
}

pub fn slash_path(path: &Path) -> String {
    let text = path.to_string_lossy().replace('\\', "/");
    text.strip_prefix("//?/").unwrap_or(&text).into()
}

fn encoded_uri(path: &Path) -> String {
    let mut encoded = String::new();
    for byte in slash_path(path).bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'/' | b':' | b'.' | b'_' | b'-' | b'~') {
            encoded.push(char::from(byte));
        } else {
            use std::fmt::Write;
            let _ = write!(encoded, "%{byte:02X}");
        }
    }
    encoded
}

fn file_uri(path: &Path) -> String {
    let encoded = encoded_uri(path);
    if encoded.starts_with("//") {
        format!("file:{encoded}")
    } else if encoded.starts_with('/') {
        format!("file://{encoded}")
    } else {
        format!("file:///{encoded}")
    }
}

pub fn sarif(diagnostics: &[Diagnostic], root: &Path) -> Value {
    let rule_ids: std::collections::BTreeSet<_> =
        diagnostics.iter().map(|item| &item.code).collect();
    let rules: Vec<_> = rule_ids.iter().map(|id| json!({"id": id})).collect();
    let results: Vec<_> = diagnostics
        .iter()
        .map(|item| {
            let location = &item.location;
            let relative = location
                .path
                .strip_prefix(root)
                .ok()
                .or_else(|| (!location.path.is_absolute()).then_some(location.path.as_path()));
            let artifact = match relative {
                Some(path) => json!({"uri": encoded_uri(path), "uriBaseId": "%SRCROOT%"}),
                None => json!({"uri": file_uri(&location.path)}),
            };
            json!({
                "ruleId": item.code,
                "level": item.severity,
                "message": {"text": item.message},
                "properties": {"subject": item.subject},
                "locations": [{"physicalLocation": {
                    "artifactLocation": artifact,
                    "region": {"startLine": location.line, "startColumn": location.column}
                }}]
            })
        })
        .collect();
    json!({
        "$schema": "https://json.schemastore.org/sarif-2.1.0.json",
        "version": "2.1.0",
        "runs": [{
            "tool": {"driver": {"name": "fortran-check", "version": env!("CARGO_PKG_VERSION"), "rules": rules}},
            "originalUriBaseIds": {"%SRCROOT%": {"uri": format!("{}/", file_uri(root).trim_end_matches('/'))}},
            "results": results
        }]
    })
}
