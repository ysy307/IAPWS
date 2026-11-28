use std::collections::BTreeSet;

use fortran_check::config::{Config, External};
use fortran_check::diagnostic::Severity;
use fortran_check::numerical;
use fortran_check::semantic::Project;
use fortran_check::source::{Source, Unit};

fn report(source: &str, config: &Config) -> BTreeSet<(String, Severity, String)> {
    let units = [Unit::parse(Source::new("numerical.F90", source)).unwrap()];
    assert!(
        units[0].diagnostics().is_empty(),
        "{:?}",
        units[0].diagnostics()
    );
    let project = Project::new(&units, config);
    numerical::check(&project, config)
        .into_iter()
        .map(|item| {
            let marker = source.lines().nth(item.location.line - 1).unwrap();
            let marker = marker
                .split_once('!')
                .map(|(_, label)| label.trim())
                .unwrap_or_default();
            (item.code, item.severity, marker.to_owned())
        })
        .collect()
}

fn expected(items: &[(&str, Severity, &str)]) -> BTreeSet<(String, Severity, String)> {
    items
        .iter()
        .map(|(code, severity, marker)| (code.to_string(), *severity, marker.to_string()))
        .collect()
}

#[test]
fn kernel_purity_and_internal_io() {
    let source = r#"module model
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none(type, external)
    private
contains
    subroutine compute_bad() ! missing purity
        implicit none(type, external)
        print *, 'stdout' ! print
        write (*, '(a)') 'stdout' ! write
    end subroutine
    pure subroutine compute_text(n)
        implicit none(type, external)
        integer(int32), intent(inout) :: n
        character(len=32) :: buffer
        write (unit=buffer, fmt='(i0)') n
        read (buffer, *) n
        inquire (iolength=n) buffer
    end subroutine
    impure elemental subroutine compute_impure(n) ! impure elemental
        implicit none(type, external)
        integer(int32), intent(inout) :: n
        n = n + 1_int32
    end subroutine
    elemental subroutine compute_elemental(n)
        implicit none(type, external)
        integer(int32), intent(inout) :: n
        n = n + 1_int32
    end subroutine
    subroutine unselected()
        implicit none(type, external)
        print *, 'allowed outside kernel'
    end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.kernels = vec!["model::compute_*".into()];
    assert_eq!(
        report(source, &config),
        expected(&[
            ("NUM001", Severity::Error, "missing purity"),
            ("NUM001", Severity::Error, "impure elemental"),
            ("NUM002", Severity::Error, "print"),
            ("NUM002", Severity::Error, "write"),
        ])
    );
}

#[test]
fn separate_body_inherits_purity_and_selection() {
    let source = r#"module model
    implicit none(type, external)
    private
    interface
        pure module subroutine compute()
            implicit none(type, external)
        end subroutine
    end interface
end module
submodule (model) model_impl
    implicit none(type, external)
contains
    module procedure compute
        print *, 'stdout' ! separate body
    end procedure
end submodule
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.kernels = vec!["model::compute".into()];
    assert_eq!(
        report(source, &config),
        expected(&[("NUM002", Severity::Error, "separate body"),])
    );
}

#[test]
fn selected_policy_warns_on_callback_purity_and_still_checks_kernel_io() {
    let source = r#"module app_ftcms_kernels
    implicit none(type, external)
    private
contains
    subroutine f0_probe() bind(C) ! callback purity
        implicit none(type, external)
        print *, 'stdout' ! callback I/O
    end subroutine
end module
module constitutive_laws
    implicit none(type, external)
    private
contains
    subroutine evaluate_probe() ! calculation purity
        implicit none(type, external)
    end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.kernels = vec![
        "app_ftcms_kernels*::[fg][0-3]_*".into(),
        "constitutive_laws::evaluate_*".into(),
    ];
    config.numerical.purity_warnings = vec!["app_ftcms_kernels*::[fg][0-3]_*".into()];
    assert_eq!(
        report(source, &config),
        expected(&[
            ("NUM001", Severity::Warning, "callback purity"),
            ("NUM002", Severity::Error, "callback I/O"),
            ("NUM001", Severity::Error, "calculation purity"),
        ])
    );
}

#[test]
fn known_real_precision_follows_aliases_and_explicit_casts() {
    let source = r#"module precision
    use, intrinsic :: iso_fortran_env, only: promoted=>real64, base=>real32, int32
    implicit none(type, external)
    private
    public :: high, low
    integer(int32), parameter :: high=promoted, low=base
end module
module model
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_double, c_float
    use precision, only: wp=>high, sp=>low
    implicit none(type, external)
    private
    integer(int32), parameter :: hp=kind(1.0_wp)
    real(sp), parameter :: initial=1.0_hp ! initialization
    type :: state_type
        private
        real(wp) :: high
    end type
contains
    subroutine update(state)
        implicit none(type, external)
        type(state_type), intent(in) :: state
        real(wp) :: a, x(2)
        real(sp) :: b
        real(c_float) :: single
        real(c_double) :: double
        complex(wp) :: z
        logical :: k
        b = (state%high) ! component assignment
        b = -high_result() ! function result
        a = 1.0_hp + 2.0_sp ! mixed arithmetic
        k = a < b ! mixed comparison
        b = real(a, kind=sp)
        b = real(a, sp)
        b = real(a, kind=wp) ! cast remains wide
        b = abs(z) ! intrinsic result
        b = aimag(z) ! imaginary component
        b = z ! complex to real
        a = z ! same-width complex to real
        a = real(z)
        b = 1.0_sp ** 2_int32
        x = [real(wp) :: b, a]
        single = double ! C precision
        a = double
    end subroutine
    pure real(wp) function high_result()
        implicit none(type, external)
        high_result = 1.0_wp
    end function
end module
"#;
    let config = Config::load("generic", None).unwrap();
    assert_eq!(
        report(source, &config),
        expected(&[
            ("NUM003", Severity::Warning, "mixed arithmetic"),
            ("NUM003", Severity::Warning, "mixed comparison"),
            ("NUM004", Severity::Error, "initialization"),
            ("NUM004", Severity::Error, "component assignment"),
            ("NUM004", Severity::Error, "function result"),
            ("NUM004", Severity::Error, "cast remains wide"),
            ("NUM004", Severity::Error, "intrinsic result"),
            ("NUM004", Severity::Error, "imaginary component"),
            ("NUM004", Severity::Error, "complex to real"),
            ("NUM004", Severity::Error, "same-width complex to real"),
            ("NUM004", Severity::Error, "C precision"),
        ])
    );
}

#[test]
fn petsc_reexported_kind_after_cpp() {
    let source = r#"module model
    use, intrinsic :: iso_fortran_env, only: real32
    use petscsys, only: petsc_double => c_double
    implicit none(type, external)
    private
contains
    subroutine update()
        implicit none(type, external)
        real(real32) :: low
        real(kind=petsc_double) :: high
        low = high ! PETSc double precision
    end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.externals.insert(
        "petscsys".into(),
        External {
            symbols: vec!["c_double".into()],
            ..External::default()
        },
    );
    config
        .numerical
        .kinds
        .insert("petscsys::c_double".into(), 64);
    assert_eq!(
        report(source, &config),
        expected(&[("NUM004", Severity::Error, "PETSc double precision"),])
    );
}

#[test]
fn unordered_kinds_and_shadowed_intrinsic() {
    let source = r#"module model
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32
    use, intrinsic :: iso_c_binding, only: c_double, c_long_double
    implicit none(type, external)
    private
    integer(int32), parameter :: local_kind = 12
contains
    subroutine update()
        implicit none(type, external)
        real(real32) :: low
        real(real64) :: high
        real(local_kind) :: processor_specific
        real(c_long_double) :: c_extended
        real(c_double) :: c_double_value
        low = real(high, real32) ! user resolves to higher kind
        low = processor_specific
        c_double_value = c_extended
        high = high + processor_specific
    end subroutine
    pure real(real64) function real(x, opt) result(y)
        implicit none(type, external)
        real(real64), intent(in) :: x
        integer(int32), intent(in) :: opt
        y = x
    end function
end module
"#;
    let config = Config::load("generic", None).unwrap();
    assert_eq!(
        report(source, &config),
        expected(&[("NUM004", Severity::Error, "user resolves to higher kind"),])
    );
}

#[test]
fn hot_loop_dispatch_and_allocation() {
    let source = r#"module model
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none(type, external)
    private
    abstract interface
        subroutine action()
            implicit none(type, external)
        end subroutine
    end interface
    type :: state_type
        procedure(action), pointer, nopass :: pointer_member
    contains
        procedure, public :: update => apply
        procedure, public, non_overridable :: fixed => apply
        generic, public :: run => update
    end type
contains
    subroutine step(self, exact, callback)
        implicit none(type, external)
        class(state_type), intent(inout) :: self
        type(state_type), intent(inout) :: exact
        procedure(action), pointer :: callback
        integer(int32), allocatable :: a(:)
        integer(int32) :: i
        character(:), allocatable :: text
        character(len=10) :: fixed_text
        call self%update()
        allocate(a(2))
        do i = 1, 2
            call self%update() ! polymorphic dispatch
            call self%run() ! generic dispatch
            call self%fixed()
            call exact%update()
            call exact%pointer_member() ! component pointer
            call callback() ! local pointer
            allocate(a(i)) ! allocation
            a = [i] ! possible reallocation
            a = produce(i) ! allocatable result
            a(:) = [i]
            a = a + i
            a = 0_int32
            text = trim(fixed_text) ! deferred length
            text = text
        end do
    end subroutine
    subroutine apply(self)
        implicit none(type, external)
        class(state_type), intent(inout) :: self
    end subroutine
    function produce(n) result(values)
        implicit none(type, external)
        integer(int32), intent(in) :: n
        integer(int32), allocatable :: values(:)
        allocate(values(n))
    end function
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.hot_loops = vec!["model::step".into()];
    assert_eq!(
        report(source, &config),
        expected(&[
            ("NUM005", Severity::Warning, "polymorphic dispatch"),
            ("NUM005", Severity::Warning, "generic dispatch"),
            ("NUM005", Severity::Warning, "component pointer"),
            ("NUM005", Severity::Warning, "local pointer"),
            ("NUM006", Severity::Warning, "allocation"),
            ("NUM006", Severity::Warning, "possible reallocation"),
            ("NUM006", Severity::Warning, "allocatable result"),
            ("NUM006", Severity::Warning, "deferred length"),
        ])
    );
}

#[test]
fn select_type_preserves_polymorphism_only_in_class_branch() {
    let source = r#"module model
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none(type, external)
    private
    type :: state_type
    contains
        procedure, public :: update
    end type
contains
    subroutine step(x)
        implicit none(type, external)
        class(*), intent(inout) :: x
        integer(int32) :: i
        select type (obj => x)
        class is (state_type)
            do i = 1, 3
                call obj%update() ! class branch
            end do
        type is (state_type)
            do i = 1, 3
                call obj%update()
            end do
        end select
    end subroutine
    subroutine update(self)
        implicit none(type, external)
        class(state_type), intent(inout) :: self
    end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.hot_loops = vec!["model::step".into()];
    assert_eq!(
        report(source, &config),
        expected(&[("NUM005", Severity::Warning, "class branch"),])
    );
}

#[test]
fn custom_io_and_loop_syntax_selectors_change_only_the_selected_checks() {
    let source = r#"module math
contains
pure subroutine advance()
    integer :: i
    integer, allocatable :: workspace(:)
    do i = 1, 2
        allocate(workspace(i)) ! default loop
        print *, i ! omitted IO node
        write (*, *) i ! retained IO node
    end do
end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.numerical.kernels = vec!["MATH::ADVANCE".into()];
    config.numerical.hot_loops = vec!["math::*".into()];
    config.numerical.hot_loop_nodes = vec!["while_statement".into()];
    config.numerical.io_statements = vec!["write_statement".into()];
    assert_eq!(
        report(source, &config),
        expected(&[("NUM002", Severity::Error, "retained IO node"),])
    );
}

#[test]
fn configured_external_widths_follow_renamed_symbols() {
    let source = r#"module math
use number_model, only: narrow => low_kind, wide => high_kind
contains
subroutine calculate()
    real(narrow) :: a
    real(wide) :: b
    a = b ! configured width
end subroutine
end module
"#;
    let mut config = Config::load("generic", None).unwrap();
    config.externals.insert(
        "number_model".into(),
        External {
            symbols: vec!["low_kind".into(), "high_kind".into()],
            ..External::default()
        },
    );
    let kinds = &mut config.numerical.kinds;
    kinds.insert("number_model::low_kind".into(), 48);
    kinds.insert("number_model::high_kind".into(), 96);
    assert_eq!(
        report(source, &config),
        expected(&[("NUM004", Severity::Error, "configured width"),])
    );
    config
        .numerical
        .kinds
        .insert("number_model::low_kind".into(), 96);
    assert!(report(source, &config).is_empty());
}
