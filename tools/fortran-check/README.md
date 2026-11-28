# FortranCheck

FortranCheck is a standalone Rust CLI for modern Fortran code and project architecture. It ships with a generic policy and an FTCMS policy. The FTCMS profile selects Fortran 2018, loads the FTCMS layer and external-module catalog, and enables Fortitude 0.10.1. The generic profile can analyze expanded free-form Fortran with the executable alone.

## Build and run

Build with Rust 1.88 or later, Cargo, and a C compiler for the generated Tree-sitter grammar. Cargo uses the checked-in dependency lockfile.

```sh
cargo install --path ./fortran-check --locked
fortran-check --version
fortran-check rules
fortran-check config --profile ftcms > fortran-check.effective.toml
```

Run these commands from the FTCMS root. Use `cargo build --manifest-path fortran-check/Cargo.toml --locked --release` to obtain the binary under `fortran-check/target/release/` instead of installing it. The executable embeds both profiles and the Fortran intrinsic-module catalog. A deployment that enables Fortitude or processes raw Fypp templates also needs the selected Fortitude or Fypp executable.

Pass the source files or directories that form one Fortran project. A relative path is resolved against `--root`. Supply all modules needed for project-wide resolution.

```sh
fortran-check check --root /path/to/project src app \
    --config rules.toml --reports-dir build/quality
```

For an existing FTCMS build, generate the Fypp files with CMake before using its compilation database. Run the command in the environment that provides the compiler, the configured libraries and the executable paths recorded in that database.

```sh
fortran-check check --root . --profile ftcms \
    --compile-db CMakeBuild/compile_commands.json --compiler-id IntelLLVM \
    --fortitude .venv/bin/fortitude --reports-dir build/quality
```

Use `--compiler-id GNU` for gfortran. For raw inputs that use CPP, add `--cpp-executable gfortran --compiler-id GNU` and repeat `--cpp-arg=-DNAME=value` or `--cpp-arg=-I/path/to/headers` for each project-specific argument. For raw `.fypp` inputs, the checker expands Fypp before CPP.

The exit codes are 0 for a completed check with no policy errors, 1 for policy errors, and 2 for incomplete analysis or tool failure. Warnings and notes appear in the output. Use `--format json` for machine-readable console output. With `--reports-dir`, the default outputs are `diagnostics.sarif`, `metrics.json`, and `dependencies.dot` for Graphviz.

## CMake integration

Include the supplied module in the consuming project. Call `fortran_check_enable` after `project(... LANGUAGES Fortran)` and before declaring compilation targets. Call `fortran_check_add_target` after declaring targets.

```cmake
include("/path/to/fortran-check/cmake/FortranCheck.cmake")
fortran_check_enable(PROFILE generic CONFIG fortran-check.toml)

add_executable(simulation src/main.F90)
fortran_check_add_target(TARGETS simulation)
```

Use `EXECUTABLE` to select a specific checker binary and `FORTITUDE` to select the Fortitude executable. `ROOT` and `REPORTS_DIR` override the project and report directories. Supply the goals that generate Fortran source files through `GENERATED_TARGETS`. A call with no `TARGETS` attaches the check to the compilation targets in the current directory and its subdirectories.

The module supports single-configuration Ninja and Makefiles generators with GNU Fortran or Intel ifx. Both Debug and Release run the same project check through normal `cmake --build`, including when a compilation target is selected directly. CMake executes the analysis after the named generation targets and before any attached compilation target. The compiler launcher requires a successful check for that build. The CMake module also applies the language and warning flags from the selected policy. The main FTCMS build uses its existing integration unless the consuming CMake project includes this module.

```sh
cmake --build build --target fortran_check
```

## Configuration and tests

See [Configuration](docs/configuration.md) for policy tables and [Rules](docs/rules.md) for diagnostic conditions and default severities. `fortran-check config --profile ftcms --config rules.toml` shows the effective policy after merging defaults, the profile and local settings.

```sh
cargo test --manifest-path fortran-check/Cargo.toml --locked
```

The CMake tests run when the relevant compiler and generator are available. Set `FORTRAN_CHECK_TEST_GFORTRAN` and `FORTRAN_CHECK_TEST_IFX` to require specific compiler paths. Set `FORTRAN_CHECK_TEST_FORTITUDE` to enable the Fortitude integration test with that executable.

## Roadmap coverage

| Phase | Implementation |
| --- | --- |
| 1 | Native language rules, Fortitude and GNU/ifx strict compilation flags |
| 2 | Encapsulation, published APIs, bindings and qualified exceptions |
| 3 | Project symbol tables, imports and renames, generics, inheritance and separate submodules |
| 4 | Inheritance, cohesion and coupling metrics and call/component graphs |
| 5 | Layer and module dependency enforcement, fan-in/out, instability and DOT |
| 6 | Pointer association, allocation, resource ownership and finalizers |
| 7 | Numerical kernels, precision, hot loops, MPI/PETSc and parallel side effects |
| 8 | Native preprocessing, source maps, build gate, caches and reports |
| 9 | Rule tests, regressions, compiler integration tests, pinned dependencies and rule reference |
| 10 | Installable standalone CLI, built-in profiles and effective policy export |
