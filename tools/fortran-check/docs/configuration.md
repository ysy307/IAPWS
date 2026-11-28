# Configuration

Use `fortran-check config --profile generic` or `--profile ftcms` to print all effective settings, including the MPI and PETSc signatures. The built-in profile is applied over the defaults. A local file passed to `--config` overlays that result. Tables merge recursively. Scalars and arrays replace the corresponding settings. Unknown keys and invalid policy values are errors. The policy format is version 1.

Paths use `/` and normally match relative to `--root`. For sources outside the root, path filters see the absolute path. Path and name patterns are case-insensitive. In paths, `*` matches within one segment and `**` can span segments. The single pattern `"*"` matches the entire value. Use qualified Fortran names such as `module::type::component` and `module::procedure::argument` as subject patterns.

## Rule selection, severity and exceptions

Every native diagnostic ID appears in `fortran-check rules`. Use a `[rules.ID]` table to set `severity` to `error`, `warning`, `note` or `off`. Optional arrays `paths` and `subjects` restrict the scope. `exclude_paths` and `exclude_subjects` exclude matched diagnostics. Empty inclusion lists accept all values.

Use `[[overrides]]` for ordered rules with `rules`, optional `paths` and `subjects`, and a `severity`. An override can select a group using a glob such as `"PRL*"`. Later matching overrides take priority among diagnostics within the original rule scope.

```toml
[rules.CBO001]
severity = "error"
paths = ["src/**"]
exclude_subjects = ["legacy_data::*"]

[[overrides]]
rules = ["NUM001"]
paths = ["src/app/**"]
subjects = ["app_ftcms_kernels*::[fg][0-3]_*"]
severity = "warning"
```

An exception requires a rule ID, an exact Fortran subject, a category, a concrete justification and optionally `paths`. Use the subject from JSON output. An `INTENT(OUT)` exception requires `significant_benefit`. Record the measured or otherwise substantial benefit and why the required argument semantics justify it.

```toml
[[exceptions]]
rule = "STD002"
subject = "io_buffer::read_values::values"
category = "significant_benefit"
reason = "Each call replaces a large allocatable buffer. Deallocation on entry is required by the caller's replacement contract and removes a measured extra copy."
paths = ["src/io/**"]
```

| Rules | Built-in categories |
| --- | --- |
| `STD002` | `significant_benefit` |
| `OOP002` | `interop`, for example an exposed component of a `bind(c)` type |
| `OOP005`, `OOP007` | `constructor`, `factory`, `callback`, `interop` |
| `OOP008` | `polymorphism`, `callback`, `interop` |
| `OOP009` | `callback`, `interop`, `method_without_object` |
| `OOP010` | `borrowed_pointer`, `interop` |
| `COH001` | `cohesion` |

Set `rules.ID.exception_categories` to customize the approved categories for a rule. Use an empty list to disallow exceptions for that rule. Use `["*"]` to accept any stated category for a rule that permits customized categories. `STD002` keeps the `significant_benefit` requirement. Parser severity settings affect the presentation of diagnostics. An incomplete parse always fails a build gate.

## Language, graph and numeric policy

| Table | Settings |
| --- | --- |
| `language` | `standard`, `line_length`, `implicit_scopes`, `require_implicit_external`, intrinsic types in `named_kinds`, and `banned_statements` |
| `limits` | `dit`, `lcom4`, `cbo`, `wmc`, and per-procedure `complexity` thresholds |
| `metrics` | `exclude_abstract` and `exclude_stateless` control LCOM4 enforcement. `external_fan_out` includes external dependencies. `complexity_nodes` selects Tree-sitter node kinds that increase complexity. |
| `layers.<name>` | `paths` assigns modules to a layer. `may_use` allows dependencies on named layers, including the same layer when listed. |
| `internals` | An array of tables with a `module` name or pattern and allowed consumer `layers`. |
| `externals.<module>` | `layers` and `paths` allow access from matching consumers. `intrinsic` identifies compiler intrinsic modules. `symbols` lists known exports. Set `compiler = true` with an empty `symbols` list to defer exports to the compiler. |
| `numerical` | `kernels`, `hot_loops`, and `purity_warnings` select qualified procedures. `io_statements` and `hot_loop_nodes` select grammar node kinds. `kinds` maps known kind constants to precision in bits. |

```toml
[limits]
dit = 3
lcom4 = 1
cbo = 5
wmc = 50
complexity = 10

[layers.physics]
paths = ["src/physics/**"]
may_use = ["math", "physics"]

[layers.math]
paths = ["src/math/**"]
may_use = ["math"]

[externals.mpi_f08]
compiler = true
layers = ["physics"]
paths = ["src/physics/mpi_adapter.F90"]

[numerical]
kernels = ["constitutive::*"]
hot_loops = ["assembly::*"]
purity_warnings = ["constitutive::recording_callback"]

[numerical.kinds]
"precision::wide_real" = 128
```

Set `language.require_implicit_external = false` to accept bare `IMPLICIT NONE`
or `IMPLICIT NONE (TYPE)` in the configured scopes. This supports compilers
without `IMPLICIT NONE (EXTERNAL)` while still checking implicit typing.

Use `limits.subjects` to set a ceiling for a specific type or procedure. Each
entry uses an exact qualified name and requires a reason. Any omitted metrics
use their global ceilings. Reports include the measured values for every
subject, and increases above a subject's ceiling still produce diagnostics.

```toml
[limits.subjects."model::state"]
cbo = 7
wmc = 53
reason = "The public state composes phase models and exposes their properties."

[limits.subjects."model::solve"]
complexity = 12
reason = "The solver uses separate guards for phase and bracket boundaries."
```

The FTCMS profile scopes the `PURE` requirement to selected numerical kernels. Residual and Jacobian callbacks that record errors for later reporting start at warning severity. The generic profile starts with empty kernel, hot-loop and layer selections.

## Resource and parallel policy

Argument selectors use a zero-based `index`. A negative index counts back from the last positional argument. The optional `keywords` array matches explicitly named actual arguments first. Use `enabled = false` to disable an existing selector or signature.

| Setting under `safety` | Purpose |
| --- | --- |
| `mpi_modules`, `petsc_modules` | Module patterns that determine API origin |
| `mpi_prefixes`, `petsc_prefixes`, `call_families` | Well-known routine naming families and exact `module::remote_name` family overrides for libraries or project wrappers |
| `resources.<acquire>` | `family`, `release`, acquisition and release handle selectors, optional status and context selectors, and whether to treat the resource as borrowed |
| `mpi_calls.<routine>` | Optional `status`, `outputs`, `rank`, `request` and `collective` signatures |
| `petsc_status`, `petsc_call_status`, `petsc_no_status` | Fallback status argument, exact per-routine status selectors and status-free PETSc calls |
| `petsc_splits.<begin>` | Matching `end` call and argument selectors for the interface's context |
| `status_handlers`, `terminators`, `finalizers`, `success_constants` | Status-consuming arguments, fatal paths, finalization calls and success values |
| `pointer_outputs` | Output pointers initialized by an API such as `c_f_pointer` |
| `self_communicators`, `request_nulls`, `collective_path_limit` | Single-rank communicator names, null-request symbols and the number of examined collective paths |
| `check_pure`, `check_do_concurrent`, `check_openmp` | Enable side-effect checks for the corresponding contexts |

A rank signature has `rank` and `communicator` selectors. A collective signature has a `communicator` selector and optional `compare_arguments`. A request signature has `action`, `handle`, and optional `count`, `completion_flag` and `collective`. The actions are `post`, `initialize`, `start`, `wait`, `test`, `cancel` and `free`. Use lower-case routine and module keys in catalogs.

```toml
[safety.call_families]
"project_petsc::buffer_create" = "petsc"
"project_petsc::buffer_destroy" = "petsc"

[safety.resources.buffer_create]
family = "petsc"
release = "buffer_destroy"
acquire_handle = { index = 1, keywords = ["buffer"] }
release_handle = { index = 0, keywords = ["buffer"] }
acquire_status = { index = 2, keywords = ["ierr"] }
release_status = { index = 1, keywords = ["ierr"] }

[safety.petsc_call_status.buffer_create]
index = 2
keywords = ["ierr"]
```

## Tools, inputs, reports and build

| Table | Settings |
| --- | --- |
| `inputs` | `extensions`, `include`, `exclude`, `discovery_exclude`, `follow_links`, and `allow_empty`. `discovery_exclude` applies when walking directories. Compiled units use scope filters on their mapped Fortran source. |
| `fortitude` | `enabled`, `executable`, pinned `version`, optional `config`, `select`, `ignore`, `preview`, `allow_comments`, `severity`, and `batch_size`. Use the FortranCheck `[rules.ID]` path filters to scope diagnostics by original filename. C121 is provided by the native checker. |
| `compiler` | `gnu_flags`, `ifx_flags`, `ifx_windows_flags` set strict options. `gnu_cpp`, `ifx_cpp`, `ifx_windows_cpp` set the native preprocess modes. A compiler array is an ordered list of individual arguments. Flags support `{standard}` and `{intel_standard}`. |
| `fypp` | `executable`, ordered `args`, and `[[fypp.overrides]]` entries with `paths` and additional `args`. Arguments support `{root}` and `{line_length}`. Retain line markers for input locations. |
| `runtime` | `jobs`, `cache`, and optional `cache_directory`. Set `--no-cache` to rerun all analysis. |
| `reports` | `include_sarif`, `include_metrics` and `include_graph` control files emitted to `--reports-dir`. |

Configure GNU and Intel warning severity through their compiler option arrays. The checker report severities govern the Rust and Fortitude diagnostics. Files with CPP include dependencies are preprocessed again to track header changes. Independent source expansion and parsing can run in parallel. Unchanged source-local Fortitude results are reused individually. The full project report is reused when the expanded sources, position maps, policy and tool stamps agree. Changes that affect the project invalidate the semantic analysis as a whole.
