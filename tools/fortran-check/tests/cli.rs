use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    time::{Duration, UNIX_EPOCH},
};

use serde_json::Value;

const VALID: &str = "program smoke\n  implicit none(type, external)\nend program smoke\n";
const INCOMPLETE_IMPLICIT: &str = "program smoke\n  implicit none\nend program smoke\n";

fn cli(args: &[&str], root: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fortran-check"))
        .args(args)
        .current_dir(root)
        .output()
        .unwrap()
}

fn status(output: &Output, expected: i32) {
    assert_eq!(
        output.status.code(),
        Some(expected),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn json(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn policy_cache_reports_and_launcher_share_the_same_result() {
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    fs::write(root.join("main.f90"), VALID).unwrap();
    fs::write(root.join("policy.toml"), "").unwrap();
    let check = [
        "check",
        "main.f90",
        "--config",
        "policy.toml",
        "--reports-dir",
        "reports",
        "--gate",
        "reports/gate.json",
        "--format",
        "json",
    ];
    let launch = [
        "launch",
        "--gate",
        "reports/gate.json",
        "--",
        env!("CARGO_BIN_EXE_fortran-check"),
        "--version",
    ];

    status(&cli(&launch, root), 2);
    let first = cli(&check, root);
    status(&first, 0);
    assert_eq!(json(&first)["complete"], true);
    assert_eq!(json(&first)["cached"], false);
    assert_eq!(json(&cli(&check, root))["cached"], true);
    assert!(root.join("reports/diagnostics.sarif").is_file());
    assert!(root.join("reports/metrics.json").is_file());
    assert!(root.join("reports/dependencies.dot").is_file());
    status(&cli(&launch, root), 0);

    fs::write(root.join("main.f90"), INCOMPLETE_IMPLICIT).unwrap();
    let failure = cli(&check, root);
    status(&failure, 1);
    assert_eq!(json(&failure)["diagnostics"][0]["code"], "STD001");
    assert_eq!(json(&failure)["diagnostics"][0]["location"]["line"], 2);
    status(&cli(&launch, root), 2);

    fs::write(
        root.join("policy.toml"),
        "[rules.STD001]\nseverity = 'warning'",
    )
    .unwrap();
    let advisory = cli(&check, root);
    status(&advisory, 0);
    assert_eq!(json(&advisory)["cached"], false);
    assert_eq!(json(&advisory)["diagnostics"][0]["severity"], "warning");
    status(&cli(&launch, root), 0);

    fs::write(root.join("policy.toml"), "[runtime]\nunrecognized = 1").unwrap();
    let invalid_policy = cli(&check, root);
    status(&invalid_policy, 2);
    assert_eq!(json(&invalid_policy)["diagnostics"][0]["code"], "TOOL001");
    status(&cli(&launch, root), 2);
}

#[test]
fn ignored_parse_diagnostics_cannot_publish_a_successful_gate() {
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    fs::write(root.join("main.F90"), format!("#if ON\n{VALID}#endif\n")).unwrap();
    fs::write(
        root.join("policy.toml"),
        "[rules.PAR001]\nseverity = 'off'\n[rules.PAR002]\nseverity = 'off'",
    )
    .unwrap();
    let check = cli(
        &[
            "check",
            "main.F90",
            "--config",
            "policy.toml",
            "--gate",
            "gate.json",
            "--format",
            "json",
        ],
        root,
    );
    status(&check, 2);
    assert_eq!(json(&check)["complete"], false);
    let gate: Value = serde_json::from_slice(&fs::read(root.join("gate.json")).unwrap()).unwrap();
    assert_eq!(gate["passed"], false);
}

#[test]
fn malformed_fortran_fails_with_a_parser_diagnostic() {
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    fs::write(
        root.join("main.f90"),
        "program smoke\n  implicit none(type, external)\n  call broken(\n",
    )
    .unwrap();
    let result = cli(&["check", "main.f90", "--format", "json"], root);
    status(&result, 2);
    let report = json(&result);
    assert_eq!(report["complete"], false);
    assert!(
        report["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["code"] == "PAR001")
    );
}

#[test]
fn effective_policy_roundtrips_and_lists_native_rules() {
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    let rules = cli(&["rules", "--format", "json"], root);
    status(&rules, 0);
    let actual = json(&rules);
    for id in [
        "STD002", "OOP005", "SEM004", "DIT003", "ARCH006", "MEM001", "OWN002", "NUM001", "PRL006",
        "PAR002",
    ] {
        assert!(
            actual
                .as_array()
                .unwrap()
                .iter()
                .any(|rule| rule["id"] == id),
            "{id}"
        );
    }
    let profile = cli(&["config", "--profile", "ftcms"], root);
    status(&profile, 0);
    fs::write(root.join("effective.toml"), profile.stdout).unwrap();
    let roundtrip = cli(
        &[
            "config",
            "--profile",
            "ftcms",
            "--config",
            "effective.toml",
            "--format",
            "json",
        ],
        root,
    );
    status(&roundtrip, 0);
    assert_eq!(
        json(&roundtrip)["numerical"]["purity_warnings"][0],
        "app_ftcms_kernels*::[fg][0-3]_*"
    );
    assert_eq!(json(&roundtrip)["limits"]["dit"], 3);
    let options = cli(
        &["compiler-flags", "--compiler-id", "GNU", "--format", "json"],
        root,
    );
    status(&options, 0);
    assert!(
        json(&options)
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item == "-Werror")
    );
}

#[test]
fn uses_configured_fortitude_when_available() {
    let Some(executable) = std::env::var_os("FORTRAN_CHECK_TEST_FORTITUDE").map(PathBuf::from)
    else {
        return;
    };
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    fs::write(
        root.join("bad.f90"),
        VALID.replace(
            "end program smoke",
            "  character*3 :: value\n  value = 'abc'\nend program smoke",
        ),
    )
    .unwrap();
    fs::write(
        root.join("policy.toml"),
        "[fortitude]\nselect = ['OB', 'S001']",
    )
    .unwrap();
    let check = [
        "check",
        "--fortitude",
        executable.to_str().unwrap(),
        "--config",
        "policy.toml",
        "--reports-dir",
        "reports",
        "--format",
        "json",
    ];
    let result = cli(&check, root);
    status(&result, 1);
    assert!(
        json(&result)["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|item| item["code"].as_str().unwrap().starts_with("OB"))
    );
    let cache = root.join("reports/cache/fortitude");
    let entry = fs::read_dir(&cache)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let old = UNIX_EPOCH + Duration::from_secs(60);
    fs::File::open(&entry).unwrap().set_modified(old).unwrap();
    fs::write(root.join("good.f90"), VALID.replace("smoke", "second")).unwrap();
    let incremental = cli(&check, root);
    status(&incremental, 1);
    assert_eq!(json(&incremental)["sources"], 2);
    assert_eq!(json(&incremental)["cached"], false);
    assert_eq!(fs::metadata(entry).unwrap().modified().unwrap(), old);
    assert_eq!(fs::read_dir(cache).unwrap().count(), 2);
}

#[test]
fn raw_fypp_diagnostics_follow_the_template_through_cpp() {
    let Some(fypp) = std::env::var_os("FORTRAN_CHECK_TEST_FYPP").map(PathBuf::from) else {
        return;
    };
    let compiler =
        std::env::var("FORTRAN_CHECK_TEST_GFORTRAN").unwrap_or_else(|_| "gfortran".to_owned());
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    fs::write(
        root.join("template.fypp"),
        "#:set name = 'smoke'\nprogram ${name}$\n#ifdef CHECK_BRANCH\n  implicit none\n#else\n  implicit none(type, external)\n#endif\nend program ${name}$\n",
    ).unwrap();
    fs::write(
        root.join("policy.toml"),
        format!("[fypp]\nexecutable = {:?}\n", fypp.to_string_lossy()),
    )
    .unwrap();
    let result = cli(
        &[
            "check",
            "template.fypp",
            "--config",
            "policy.toml",
            "--compiler-id",
            "GNU",
            "--cpp-executable",
            &compiler,
            "--cpp-arg=-DCHECK_BRANCH",
            "--format",
            "json",
        ],
        root,
    );
    status(&result, 1);
    let data = json(&result);
    assert_eq!(data["diagnostics"][0]["code"], "STD001");
    assert_eq!(data["diagnostics"][0]["location"]["line"], 4);
    assert!(
        data["diagnostics"][0]["location"]["path"]
            .as_str()
            .unwrap()
            .ends_with("template.fypp")
    );
}
