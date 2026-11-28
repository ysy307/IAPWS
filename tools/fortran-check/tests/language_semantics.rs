use fortran_check::config::Config;
use fortran_check::language;
use fortran_check::semantic::Project;
use fortran_check::source::{Source, Unit};
use std::collections::BTreeSet;
use std::path::Path;

fn units(sources: &[&str]) -> Vec<Unit> {
    sources
        .iter()
        .enumerate()
        .map(|(i, source)| Unit::parse(Source::new(format!("unit_{i}.f90"), *source)).unwrap())
        .collect()
}

#[test]
fn private_exports_and_renames_are_resolved_across_units() {
    let units = units(&[
        "module library
          implicit none(type, external)
          private
          integer, parameter, public :: EXPORTED = 1
          integer, parameter :: hidden = 2
        end module
        ",
        "module facade
          use Library, only: value => EXPORTED
          implicit none(type, external)
          private
          public :: value
        end module
        ",
        "program main
          use facade, only: result => value
          use library, only: hidden
          implicit none(type, external)
          print *, result
        end program
        ",
    ]);
    let project = Project::new(&units, &Config::default());
    let main = project
        .scopes
        .iter()
        .position(|scope| scope.name == "main")
        .unwrap();
    let value = project.lookup(main, "RESULT").unwrap();
    assert_eq!(project.symbols[value].qualified, "library::exported");
    assert!(project.lookup(main, "hidden").is_none());
    assert_eq!(
        project
            .diagnostics
            .iter()
            .map(|d| d.code.as_str())
            .collect::<Vec<_>>(),
        ["SEM002"]
    );
}

#[test]
fn submodule_bodies_inherit_interfaces_and_link_bound_calls() {
    let units = units(&[
        "module algorithms
          implicit none(type, external)
          private
          type, public :: model
            private
            integer :: stored
          contains
            procedure, public :: apply => calculate
          end type
          interface
            module subroutine calculate(self, x)
              class(model), intent(inout) :: self
              integer, intent(in) :: x
            end subroutine
          end interface
        end module
        ",
        "submodule(algorithms) seed
          implicit none(type, external)
          integer, parameter :: offset=2
        end submodule
        ",
        "submodule(algorithms:seed) implementation
          implicit none(type, external)
        contains
          module procedure calculate
            self%stored=x+offset
          end procedure
        end submodule
        ",
    ]);
    let project = Project::new(&units, &Config::default());
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let binding = project.member(project.types[0], "apply").unwrap();
    let body = project
        .implementation(project.binding_target(binding).unwrap())
        .unwrap();
    assert_eq!(
        project.scopes[body].qualified,
        "algorithms:implementation::calculate"
    );
    assert!(project.references.iter().any(|entry| {
        entry.scope == body
            && project.symbols[entry.symbol].name == "stored"
            && entry
                .root_receiver
                .is_some_and(|id| project.symbols[id].name == "self")
    }));
}

#[test]
fn all_language_and_oop_rules_have_precise_violations() {
    let units = units(&["module bad
      ! No private default or explicit implicit rules.
      implicit none
      integer :: global = 0
      type :: state
        integer, pointer :: values(:)
      contains
        procedure :: run
        procedure, nopass, public :: cleanup
      end type
      interface factory
        module procedure absent
      end interface
    contains
      subroutine run(self, result)
        implicit none(type, external)
        class(state) :: self
        integer, intent(out) :: result
        external :: legacy
        go to 10
        10 continue
      end subroutine
      subroutine cleanup()
        implicit none(type, external)
      end subroutine
      subroutine callback(value)
        implicit none(type, external)
        class(*) :: value
      end subroutine
    end module
    "]);
    let config = Config::default();
    let project = Project::new(&units, &config);
    let lint = language::check(&project, &config);
    let codes: BTreeSet<_> = lint.iter().map(|item| item.code.as_str()).collect();
    for code in [
        "STD001", "STD002", "STD003", "STD004", "STD005", "STD006", "OOP001", "OOP002", "OOP003",
        "OOP004", "OOP005", "OOP006", "OOP007", "OOP008", "OOP009", "OOP010", "OOP011",
    ] {
        assert!(
            codes.contains(code),
            "Missing {code}. Reported {codes:?}, semantic {:?}",
            project.diagnostics
        );
    }
    assert!(
        lint.iter()
            .any(|d| d.code == "STD002" && d.subject == "bad::run::result")
    );
}

#[test]
fn abstract_interface_and_inherited_generic_preserve_private_components() {
    let units = units(&[
        "module factory
          use, intrinsic :: iso_fortran_env, only: int32
          implicit none(type, external)
          private
          type, abstract :: model
          contains
            procedure(operation), public, deferred :: execute
          end type
          abstract interface
            subroutine operation(self, value)
              import, only: model, int32
              implicit none(type, external)
              class(model), intent(inout) :: self
              integer(int32), intent(in) :: value
            end subroutine
          end interface
          type, extends(model), public :: instance
            private
            integer(int32) :: stored
          contains
            procedure, public :: execute => run
            generic, public :: update => execute
          end type
        contains
          subroutine run(self, value)
            implicit none(type, external)
            class(instance), intent(inout) :: self
            integer(int32), intent(in) :: value
            self%stored = value
          end subroutine
        end module
        ",
        "program driver
          use, intrinsic :: iso_fortran_env, only: int32
          use factory, only: instance
          implicit none(type, external)
          type(instance), save :: object
          call object%update(12_int32)
        end program
        ",
    ]);
    let config = Config::load("generic", None).unwrap();
    let project = Project::new(&units, &config);
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    assert!(language::check(&project, &config).is_empty());
    let (base, instance) = (project.types[0], project.types[1]);
    assert_eq!(project.resolve_parent(instance), Some(base));
    let binding = project
        .binding_target(project.member(base, "execute").unwrap())
        .unwrap();
    assert_eq!(project.symbols[binding].kind, "abstract_interface");
    assert!(project.references.iter().any(|reference| {
        reference.call
            && project.symbols[reference.symbol].qualified == "factory::instance::execute"
    }));
}

#[test]
fn intent_exceptions_are_exact_and_consistent_across_separate_procedures() {
    let config = Config::from_toml(
        "
        [[exceptions]]
        rule = 'STD002'
        subject = 'warehouse::fetch::left'
        category = 'significant_benefit'
        reason = 'Copies would allocate an equally sized temporary result buffer.'
        [externals.iso_fortran_env]
        intrinsic = true
        symbols = ['int32']
    ",
    )
    .unwrap();
    let units = units(&[
        "MoDuLe warehouse
          use, intrinsic :: iso_fortran_env, only: ik => int32
          implicit none(type, &
          & external)
          private
          interface
            module subroutine fetch(left, right)
              implicit none(type, external)
              integer(ik), intent (out) :: left, right
            end subroutine
          end interface
        end module
        ",
        "submodule(warehouse) implementation
          implicit none(type, external)
        contains
          module subroutine fetch(left, right)
            implicit none(type, external)
            integer(ik), intent & ! Continuation between attribute and selector.
              & (OuT) :: left, right
            left = 1_ik
            right = 2_ik
          end subroutine
        end submodule
        ",
    ]);
    let project = Project::new(&units, &config);
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let raw = language::check(&project, &config);
    assert_eq!(
        raw.iter()
            .filter(|finding| finding.code == "STD002")
            .count(),
        4
    );
    let actual: Vec<_> = raw
        .into_iter()
        .filter_map(|item| config.filter(item, Path::new(".")))
        .collect();
    assert_eq!(actual.len(), 2, "{actual:?}");
    assert!(
        actual
            .iter()
            .all(|item| item.code == "STD002" && item.subject == "warehouse::fetch::right")
    );
}

#[test]
fn selectors_for_location_intrinsics_exclude_boz_data_and_user_calls() {
    let units = units(&[
        r#"program locations
          use, intrinsic :: iso_fortran_env, only: int32
          implicit none(type, external)
          integer(int32), parameter :: values(*) = [1_int32, 2_int32]
          logical, parameter :: include_all = .true.
          print *, findloc(values, 1_int32, 1_int32, include_all, 4)
          print *, FINDLOC(values, 1_int32, include_all, 4)
          print *, minloc(values, include_all, 4)
          print *, int(z'DEAD', kind=int32), 1_"a", .true._1
          print *, maxloc(values, include_all, int32, .true.)
        end program
        "#,
        "program user_conversion
          use, intrinsic :: iso_fortran_env, only: ik => int32
          implicit none(type, external)
          print *, int(1_ik, 4)
        contains
          integer(ik) function int(first, second)
            implicit none(type, external)
            integer(ik), intent(in) :: first, second
            int = first + second
          end function
        end program
        ",
    ]);
    let config = Config::load("generic", None).unwrap();
    let project = Project::new(&units, &config);
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let findings = language::check(&project, &config);
    let lines: Vec<_> = findings
        .iter()
        .map(|item| (item.code.as_str(), item.location.line))
        .collect();
    assert_eq!(
        lines,
        [
            ("STD004", 6),
            ("STD004", 7),
            ("STD004", 8),
            ("STD004", 9),
            ("STD004", 9)
        ]
    );
}

#[test]
fn unresolved_duplicate_inaccessible_and_ambiguous_symbols_are_diagnosed() {
    let units = units(&[
        "module a
          private
          integer, parameter, public :: shared = 1
          type, public :: record
            private
            integer :: secret
          end type
          interface
            module subroutine absent()
            end subroutine
          end interface
          interface unspecified
            module procedure wrong
          end interface
        end module
        ",
        "module b
          integer, parameter :: shared = 2
          integer, parameter :: duplicate = 2, duplicate = 3
        end module
        ",
        "program demo
          use a, only: record, shared
          use b ! C121 is independent of comment-based allowances.
          type(record) :: instance
          print *, instance%secret, shared, undeclared
        end program
        ",
    ]);
    let config = Config::default();
    let project = Project::new(&units, &config);
    let codes: BTreeSet<_> = project
        .diagnostics
        .iter()
        .map(|item| item.code.as_str())
        .collect();
    assert_eq!(
        codes,
        BTreeSet::from(["SEM001", "SEM002", "SEM003", "SEM004", "SEM005"])
    );
    assert!(
        project
            .diagnostics
            .iter()
            .any(|item| item.code == "SEM003" && item.subject.ends_with("::secret"))
    );
    assert!(
        project
            .diagnostics
            .iter()
            .any(|item| item.code == "SEM003" && item.subject.ends_with("::shared"))
    );
    assert!(
        language::check(&project, &config).iter().any(|item| {
            item.code == "C121" && item.subject == "demo" && item.location.line == 3
        })
    );
}

#[test]
fn generic_overloads_use_argument_rank_for_scalars_and_sections() {
    let units = units(&["module metrics
      private
      interface length
        module procedure scalar_length, vector_length
      end interface
    contains
      integer function scalar_length(value)
        integer, intent(in) :: value
        scalar_length = value
      end function
      integer function vector_length(value)
        integer, intent(in) :: value
        dimension :: value(:)
        vector_length = size(value)
      end function
      subroutine examine()
        integer, dimension(2, 2) :: values
        print *, length(values(1, 1))
        print *, length(values(:, 1))
        print *, length([1, 2])
      end subroutine
    end module
    "]);
    let project = Project::new(&units, &Config::default());
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let calls: Vec<_> = project
        .references
        .iter()
        .filter(|item| item.call && project.scopes[item.scope].name == "examine")
        .map(|item| project.symbols[item.symbol].name.as_str())
        .collect();
    assert_eq!(calls, ["scalar_length", "vector_length", "vector_length"]);
    let examine = project
        .scopes
        .iter()
        .position(|scope| scope.name == "examine")
        .unwrap();
    let id = project.lookup(examine, "values").unwrap();
    assert_eq!(project.expression(0, project.symbols[id].node), Some(id));
}

#[test]
fn a_return_type_alone_is_an_implicit_procedure_interface() {
    let units = units(&["module contracts
      implicit none(type, external)
      private
      type :: handle
        private
        procedure(real), pointer, nopass :: callback
      end type
    contains
      subroutine run()
        implicit none(type, external)
        procedure(real), pointer :: produce
        type(logical) :: available
        available = .false.
      end subroutine
    end module
    "]);
    let config = Config::default();
    let project = Project::new(&units, &config);
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let findings: Vec<_> = language::check(&project, &config)
        .into_iter()
        .filter(|item| item.code == "STD003")
        .map(|item| item.subject)
        .collect();
    assert_eq!(
        findings,
        ["contracts::run::produce", "contracts::handle::callback"]
    );
}

#[test]
fn unused_only_imports_preserve_external_origins_after_renaming() {
    let units = units(&[
        "module wrapper
          use mpi_f08, only: local_rank => mpi_comm_rank
          private
          public :: local_rank
        end module
        ",
        "program demo
          use wrapper, only: library_rank => local_rank
        end program
        ",
    ]);
    let config = Config::load("generic", None).unwrap();
    let project = Project::new(&units, &config);
    assert!(project.diagnostics.is_empty(), "{:?}", project.diagnostics);
    let scope = project
        .scopes
        .iter()
        .position(|scope| scope.name == "demo")
        .unwrap();
    let import = project.lookup(scope, "library_rank").unwrap();
    assert_eq!(project.symbols[import].qualified, "mpi_f08::mpi_comm_rank");
    assert!(
        project
            .external_imports
            .iter()
            .any(|&(source, symbol, _)| { source == scope && symbol == import })
    );
}

#[test]
fn portable_implicit_none_still_requires_explicit_typing() {
    let units = units(&[
        "module portable\nimplicit none\nend module",
        "module typed\nimplicit none (type)\nend module",
        "module complete\nimplicit none (type, external)\nend module",
        "module external_only\nimplicit none (external)\nend module",
        "module typed_implicitly\nimplicit real(a-h)\nend module",
        "module missing\nend module",
    ]);
    let diagnostics = |config: &Config| {
        let project = Project::new(&units, config);
        language::check(&project, config)
            .into_iter()
            .filter(|item| item.code == "STD001")
            .map(|item| item.subject)
            .collect::<BTreeSet<_>>()
    };
    let config = Config::from_toml("[language]\nrequire_implicit_external = false\n").unwrap();
    assert_eq!(
        diagnostics(&config),
        ["external_only", "missing", "typed_implicitly"]
            .map(str::to_owned)
            .into_iter()
            .collect()
    );
    assert_eq!(diagnostics(&Config::default()).len(), 5);
}
