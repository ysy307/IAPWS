use std::path::Path;

use fortran_check::config::Config;
use fortran_check::diagnostic::{Diagnostic, Location, Severity};
use fortran_check::source::{Source, Unit, walk};

#[test]
fn parser_maps_expanded_source_to_its_template() {
    let mut source = Source::new(
        "generated.F90",
        "module demo\n    implicit none\nend module\n",
    );
    source.line_map = vec![
        Some(("model.fypp".into(), 20)),
        Some(("model.fypp".into(), 28)),
    ];
    let unit = Unit::parse(source).unwrap();
    let implicit = walk(unit.tree.root_node())
        .find(|node| node.kind() == "implicit_statement")
        .unwrap();
    assert!(unit.diagnostics().is_empty());
    assert_eq!(unit.location(implicit), Location::new("model.fypp", 28, 5));
}

#[test]
fn alternate_rule_severity_and_exact_exception_are_independent() {
    let config = Config::from_toml(
        r#"
[rules.STD002]
severity = "error"

[[overrides]]
rules = ["STD002"]
paths = ["src/legacy/**"]
severity = "warning"

[[exceptions]]
rule = "STD002"
subject = "model::replace_values::buffer"
category = "significant_benefit"
reason = "The old allocatable buffer is discarded on entry."
"#,
    )
    .unwrap();
    let root = Path::new(".");
    let warning = Diagnostic::new(
        "STD002",
        "intent(out)",
        Location::new("src/legacy/model.f90", 5, 1),
    )
    .with_subject("model::replace_values::count");
    assert_eq!(
        config.filter(warning, root).unwrap().severity,
        Severity::Warning
    );
    let exception = Diagnostic::new(
        "STD002",
        "intent(out)",
        Location::new("src/model.f90", 3, 1),
    )
    .with_subject("model::replace_values::buffer");
    assert!(config.filter(exception, root).is_none());
}

#[test]
fn configuration_rejects_typos_and_unjustified_intent_exceptions() {
    assert!(Config::from_toml("[numerical]\n hot_lops = []").is_err());
    assert!(Config::from_toml("[rules.STD099]\n severity = 'off'").is_err());
    let policy = r#"
[[exceptions]]
rule = "STD002"
subject = "model::replace_values::buffer"
reason = "The old allocation is discarded."
category = "interop"
"#;
    assert!(Config::from_toml(policy).is_err());
}

#[test]
fn policy_can_define_an_exact_exception_category_for_any_design_rule() {
    let policy = r#"
[rules.CC001]
exception_categories = ["generated"]

[[exceptions]]
rule = "CC001"
subject = "math::operator_remap"
category = "generated"
reason = "One generated case is needed for each supported coordinate system."
"#;
    let config = Config::from_toml(policy).unwrap();
    let message = Diagnostic::new(
        "CC001",
        "too many branches",
        Location::new("math.f90", 15, 1),
    )
    .with_subject("MATH::OPERATOR_REMAP");
    assert!(config.filter(message, Path::new(".")).is_none());
    assert!(
        Config::from_toml(&policy.replace(
            "subject = \"math::operator_remap\"",
            "subject = \"math::*\""
        ))
        .is_err()
    );
}
