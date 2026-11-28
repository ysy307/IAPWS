# Rules

FortranCheck assigns each rule a stable ID. A policy can change the severity, source scope and subject scope of each diagnosis. A conditional default uses an error for definite evidence and a warning for possible evidence or a selected advisory scope. The [configuration reference](configuration.md) describes severity overrides and documented exceptions.

## 1. Language and conventions

The FTCMS profile selects Fortran 2018. GNU Fortran and ifx provide language conformance checks, extension and procedure-interface diagnostics, and warnings treated as errors in the CMake build. Fortitude selects `E`, `OB`, `PORT021` and `S001` by default for the FTCMS profile. Consult `fortitude explain` for its additional rules.

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| STD001 | Error | A selected scope lacks `implicit none (type, external)`. A separate module procedure can inherit it from its declaration. |
| STD002 | Error | A dummy has `intent(out)` and lacks a qualified, significant-benefit exception. |
| STD003 | Error | A procedure declaration or a call relies on an implicit interface. |
| STD004 | Error | A selected intrinsic declaration, real literal, or explicit kind selector lacks a named kind constant. |
| STD005 | Error | Scope or initialization gives a variable implicit `save`. |
| STD006 | Error | A statement matches a configured prohibition. `go to` is prohibited by default. |
| C121 | Error | A `use` statement lacks an explicit `only` list. |

## 2. Object-oriented design

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| OOP001 | Error | A module's default accessibility is public. |
| OOP002 | Error | A module-defined type exposes a data component or procedure pointer. Type parameters and compiler-provided parent components are excluded. |
| OOP003 | Error | A module or submodule declares mutable data or a mutable procedure pointer. |
| OOP004 | Error | A specific or generic type binding omits an explicit access attribute. |
| OOP005 | Error | A conventional module procedure or external procedure is publicly available without an approved exception. |
| OOP006 | Error | The implementation of a type-bound or final procedure is exported directly. |
| OOP007 | Error | A public module generic lacks an approved exception. |
| OOP008 | Error | `class(*)` lacks an approved exception. |
| OOP009 | Error | `nopass` lacks an approved exception. |
| OOP010 | Error | A data pointer component lacks an approved exception. |
| OOP011 | Error | A public generic contains an invalid, missing, or repeated specific. |

## 3. Project semantics

The resolver handles mixed case, `USE` renames and re-exports, host association, procedure interfaces, derived types, inheritance, generics and separate implementations. The compiler checks exact interface and overload compatibility. The configured export catalog covers declared external symbols. Compiler-managed external modules rely on the compiler for their installed exports.

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| SEM001 | Error | A declaration duplicates a symbol in its namespace. |
| SEM002 | Error | A symbol or imported name cannot be resolved. |
| SEM003 | Error | A reference is inaccessible or its candidates are ambiguous. |
| SEM004 | Error | A submodule or separate module procedure cannot be linked to its declaration or ancestor. |
| SEM005 | Error | A binding, generic specific or explicit procedure interface cannot be resolved. |

## 4. Inheritance, cohesion and coupling

| ID or metric | Default | Computation |
| --- | --- | --- |
| DIT001 | Error | The resolved extension graph contains a cycle. |
| DIT002 | Error | The base type or its depth cannot be resolved. |
| DIT003 | Error above 3 | DIT counts extension edges from the root type to the current type. |
| NOC | Metric | Count the type's direct children. |
| COH001 | Error above 1 | LCOM4 counts connected components among concrete local bindings. Two methods connect through shared passed-object state or a method call. Abstract and stateless types are exempt from this threshold by default. |
| CBO001 | Warning above 5 | Count distinct other derived types linked through inheritance, declarations and method references. |
| WMC001 | Warning above 50 | Sum cyclomatic complexity for the type's concrete methods. |
| CC001 | Warning above 10 | Start at 1 and add one for each configured control-flow node in a procedure. |
| RFC | Metric | Count the union of analyzed methods and distinct direct call targets. |

The JSON report records method calls and component accesses along with the aggregate metrics.

## 5. Architecture

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| ARCH001 | Error | A module or program has no configured layer or matches multiple layers when layer policy is active. |
| ARCH002 | Error | A project module has multiple definitions. |
| ARCH003 | Error | An imported module is absent from the project and the external catalog. |
| ARCH004 | Error | A module or program depends on a forbidden project layer. |
| ARCH005 | Error | An internal module is accessed from outside its allowed layers. |
| ARCH006 | Error | An external module is accessed outside its permitted layers or source paths. |
| ARCH007 | Error | A cycle connects project modules. |

The module graph uses importer-to-dependency edges. Fan-in counts distinct project dependents \(C_a\). Fan-out counts distinct dependencies \(C_e\), with external imports included by default. The instability metric is \(C_e / (C_a + C_e)\), with a value of zero for an isolated module. Use `dependencies.dot` to render the graph.

## 6. Memory and ownership

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| MEM001 | Warning | A locally owned pointer has usage that can be represented by an allocatable variable. |
| MEM002 | Conditional | A pointer is used with a potentially undefined association status. |
| MEM003 | Conditional | A pointer can reference null or released storage, or is deallocated while aliasing an unsuitable target. |
| MEM004 | Warning | Whole allocatable assignment can reallocate storage. |
| MEM005 | Warning | A strided actual array can require packing for a resolved dummy argument. |
| OWN001 | Warning | A configured acquired resource can reach a normal exit without its matching release. |
| OWN002 | Conditional | A release lacks a matching acquisition or uses a different access context. |
| OWN003 | Warning | Intrinsic assignment can copy a known owned value or acquired handle. |
| OWN004 | Conditional | A finalizer is invalid or an owned component has no recognized cleanup. |
| OWN005 | Warning | An acquisition or assignment can replace an active resource. |

Resource tracking follows recognized calls and local control flow. Register project wrappers and library-specific acquisition pairs in `safety`. Calls with opaque side effects and aliased storage outside the analyzed procedure can require review with the compiler or runtime tooling.

## 7. Numerical and parallel code

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| NUM001 | Conditional | A selected numerical kernel lacks `pure`. Configured callbacks can use the advisory default. |
| NUM002 | Error | A selected numerical kernel uses external I/O. |
| NUM003 | Warning | An expression or assignment mixes reals of distinct known precisions. |
| NUM004 | Error | An assignment or resolved argument conversion implicitly narrows a known precision. |
| NUM005 | Warning | A selected hot loop contains a call with dynamic type dispatch. |
| NUM006 | Warning | A selected hot loop performs allocation. |
| PRL001 | Conditional | A recognized MPI call lacks expected status handling or loses its status before handling. |
| PRL002 | Conditional | A recognized PETSc call lacks expected status handling or loses its status before handling. |
| PRL003 | Warning | A path leaves an MPI operation or PETSc Begin/End sequence unfinished. |
| PRL004 | Conditional | A request is reused, released or operated on in a state inconsistent with its communication lifecycle. |
| PRL005 | Warning | Rank-dependent paths can issue different collectives or arguments on a tracked communicator. |
| PRL006 | Conditional | A recognized pure or parallel region writes nonlocal or shared static storage. |

Numerical precision checks use known kinds and the configured precision catalog. The compiler determines implementation-specific kinds. Collective path analysis is bounded and follows visible control flow and config-defined signatures. Cross-rank behavior also depends on runtime data and communication performed through libraries.

## 8. Input and execution failures

| ID | Default | Diagnostic condition |
| --- | --- | --- |
| PAR001 | Error | The pinned Fortran grammar reports incomplete or unparseable syntax. |
| PAR002 | Error | Unexpanded preprocessing directives or includes remain in the parsed source. |
| TOOL001 | Fatal | An execution, configuration, preprocessing, input or reporting failure prevents a complete result. |

An incomplete input fails the compiler gate as well as the standalone check. `TOOL001` is an execution status and has no rule severity override. The syntax engine targets free-form Fortran.
