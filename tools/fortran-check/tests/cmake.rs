use std::{
    env, fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use serde_json::Value;

fn command(variable: &str, default: &str) -> Option<PathBuf> {
    let explicit = env::var_os(variable);
    let path = explicit
        .clone()
        .map_or_else(|| PathBuf::from(default), PathBuf::from);
    let available = Command::new(&path)
        .arg("--version")
        .output()
        .is_ok_and(|output| output.status.success());
    assert!(
        explicit.is_none() || available,
        "{variable} is set but {} is unavailable",
        path.display()
    );
    available.then_some(path)
}

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().replace('\\', "/")
}

fn require(output: Output, pass: bool) -> String {
    let out = String::from_utf8_lossy(&output.stdout);
    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        output.status.success(),
        pass,
        "stdout:\n{out}\nstderr:\n{err}"
    );
    format!("{out}\n{err}")
}

fn scenario(variable: &str, default: &str, build_type: &str, generator: &str) {
    let Some(cmake) = command("FORTRAN_CHECK_TEST_CMAKE", "cmake") else {
        return;
    };
    if generator == "Ninja" && command("FORTRAN_CHECK_TEST_NINJA", "ninja").is_none() {
        return;
    }
    if generator == "Unix Makefiles" && command("FORTRAN_CHECK_TEST_MAKE", "make").is_none() {
        return;
    }
    let Some(compiler) = command(variable, default) else {
        return;
    };
    let app = text(env!("CARGO_BIN_EXE_fortran-check"));
    let integration = text(Path::new(env!("CARGO_MANIFEST_DIR")).join("cmake/FortranCheck.cmake"));
    let base = tempfile::tempdir().unwrap();
    let root = base.path();
    let build = root.join("build");
    let cmakelists = format!(
        r#"cmake_minimum_required(VERSION 3.20)
project(gate LANGUAGES Fortran)
include("{integration}")
fortran_check_enable(EXECUTABLE "{app}" CONFIG policy.toml)
add_custom_command(
    OUTPUT "${{CMAKE_CURRENT_BINARY_DIR}}/unit.F90"
    COMMAND "${{CMAKE_COMMAND}}" -E copy
        "${{CMAKE_CURRENT_SOURCE_DIR}}/template.F90" "${{CMAKE_CURRENT_BINARY_DIR}}/unit.F90"
    DEPENDS "${{CMAKE_CURRENT_SOURCE_DIR}}/template.F90"
    VERBATIM
)
add_custom_target(generated DEPENDS "${{CMAKE_CURRENT_BINARY_DIR}}/unit.F90")
add_executable(smoke "${{CMAKE_CURRENT_BINARY_DIR}}/unit.F90")
target_include_directories(smoke PRIVATE "${{CMAKE_CURRENT_SOURCE_DIR}}")
fortran_check_add_target(GENERATED_TARGETS generated)
"#,
    );
    fs::write(root.join("CMakeLists.txt"), cmakelists).unwrap();
    fs::write(root.join("policy.toml"), "").unwrap();
    let source = root.join("template.F90");
    fs::write(
        &source,
        format!("#include \"switch.inc\"\n#line 2 \"{}\"\nprogram smoke\n#ifdef SAFE\n  implicit none(type, external)\n#else\n  implicit none\n#endif\nend program smoke\n", text(&source)),
    ).unwrap();
    fs::write(root.join("switch.inc"), "#define SAFE\n").unwrap();
    let run = |args: &[&str]| {
        Command::new(&cmake)
            .args(args)
            .current_dir(root)
            .output()
            .unwrap()
    };
    let compiler = format!("-DCMAKE_Fortran_COMPILER={}", text(compiler));
    require(
        run(&[
            "-S",
            ".",
            "-B",
            "build",
            "-G",
            generator,
            &compiler,
            &format!("-DCMAKE_BUILD_TYPE={build_type}"),
        ]),
        true,
    );
    let success = require(run(&["--build", "build", "--parallel", "4"]), true);
    assert!(
        success.contains("FortranCheck: 1 sources, 0 errors, 0 warnings"),
        "{success}"
    );
    let compiled = build.join(if cfg!(windows) { "smoke.exe" } else { "smoke" });
    assert!(compiled.is_file());
    let clean = require(run(&["--build", "build", "--parallel", "4"]), true);
    assert!(clean.contains("(cached)"), "{clean}");

    fs::write(root.join("switch.inc"), "! SAFE is undefined\n").unwrap();
    let failure = require(
        run(&["--build", "build", "--target", "smoke", "--parallel", "4"]),
        false,
    );
    assert!(
        failure
            .replace('\\', "/")
            .contains(&format!("{}:6:", text(&source)))
            && failure.contains("[STD001]"),
        "{failure}"
    );
    let gate = build.join("fortran-check/gate.json");
    let read_gate = || serde_json::from_slice::<Value>(&fs::read(&gate).unwrap()).unwrap();
    assert_eq!(read_gate()["passed"], false);

    fs::write(
        root.join("policy.toml"),
        "[rules.STD001]\nseverity = 'warning'\n",
    )
    .unwrap();
    let warning = require(run(&["--build", "build", "--parallel", "4"]), true);
    assert!(warning.contains("0 errors, 1 warnings"), "{warning}");
    assert_eq!(read_gate()["passed"], true);
    let db = fs::read_to_string(build.join("compile_commands.json")).unwrap();
    if default == "ifx" {
        assert!(db.contains("f18") && db.contains("stderrors"), "{db}");
    } else {
        assert!(db.contains("-std=f2018") && db.contains("-Werror"), "{db}");
    }

    fs::write(
        root.join("policy.toml"),
        "[rules.STD004]\nseverity = 'off'\n",
    )
    .unwrap();
    fs::write(
        source,
        "program smoke\n  implicit none(type, external)\n  real*8, save :: legacy\n  legacy = 1\n  print *, legacy\nend program smoke\n",
    )
    .unwrap();
    let rejected = require(run(&["--build", "build", "--parallel", "4"]), false);
    assert!(
        rejected.contains("FortranCheck: 1 sources, 0 errors, 0 warnings")
            && rejected.contains("real*8"),
        "{rejected}"
    );
    assert_eq!(read_gate()["passed"], true);
}

#[test]
fn gfortran_debug_ninja() {
    scenario("FORTRAN_CHECK_TEST_GFORTRAN", "gfortran", "Debug", "Ninja");
}

#[test]
fn gfortran_release_make() {
    if !cfg!(windows) {
        scenario(
            "FORTRAN_CHECK_TEST_GFORTRAN",
            "gfortran",
            "Release",
            "Unix Makefiles",
        );
    }
}

#[test]
fn ifx_debug_ninja() {
    scenario("FORTRAN_CHECK_TEST_IFX", "ifx", "Debug", "Ninja");
}

#[test]
fn ifx_release_ninja() {
    scenario("FORTRAN_CHECK_TEST_IFX", "ifx", "Release", "Ninja");
}
