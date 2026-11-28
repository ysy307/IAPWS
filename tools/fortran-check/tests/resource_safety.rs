use fortran_check::{
    config::Config,
    diagnostic::{Diagnostic, Severity},
    memory, ownership, parallel,
    semantic::Project,
    source::{Source, Unit},
};

fn diagnose(input: &str, config: &Config) -> Vec<Diagnostic> {
    let units = [Unit::parse(Source::new("safety.f90", input)).unwrap()];
    assert!(
        units[0].diagnostics().is_empty(),
        "{:?}",
        units[0].diagnostics()
    );
    let project = Project::new(&units, config);
    let mut diagnostics = memory::check(&project, config);
    diagnostics.extend(ownership::check(&project, config));
    diagnostics.extend(parallel::check(&project, config));
    diagnostics
}

fn has(diagnostics: &[Diagnostic], code: &str) -> bool {
    diagnostics.iter().any(|diagnostic| diagnostic.code == code)
}

#[test]
fn argument_keyword_does_not_read_the_pointer_with_the_same_name() {
    let source = r#"
program pointers
  implicit none (type, external)
  integer, pointer :: ref
  integer, target :: value
  call initialize(ref=ref)
  print *, ref
  nullify(ref)
  call consume(ref=ref)
contains
  subroutine initialize(ref)
    implicit none (type, external)
    integer, pointer, intent(inout) :: ref
    ref => value
  end subroutine

  subroutine consume(ref)
    implicit none (type, external)
    integer, intent(in) :: ref
    print *, ref
  end subroutine
end program
"#;
    let result = diagnose(source, &Config::default());
    assert!(!has(&result, "MEM002"), "{result:?}");
    let invalid: Vec<_> = result
        .iter()
        .filter(|issue| issue.code == "MEM003")
        .collect();
    assert_eq!(invalid.len(), 1, "{result:?}");
    assert_eq!(
        (invalid[0].location.line, invalid[0].location.column),
        (9, 20)
    );
}

#[test]
fn pointer_state_distinguishes_allocation_and_dangling_aliases() {
    let input = r#"
program pointers
  implicit none (type, external)
  integer, pointer :: ptr, alias, owned
  integer :: result
  result = ptr
  allocate(ptr)
  alias => ptr
  nullify(ptr)
  result = alias
  deallocate(alias)
  deallocate(alias)
  allocate(owned)
  deallocate(owned)
end program
"#;
    let diagnostics = diagnose(input, &Config::default());
    assert!(has(&diagnostics, "MEM001"), "{diagnostics:?}");
    assert_eq!(
        diagnostics
            .iter()
            .filter(|item| item.code == "MEM002")
            .count(),
        1,
        "{diagnostics:?}"
    );
    assert!(has(&diagnostics, "MEM003"), "{diagnostics:?}");
}

#[test]
fn status_terminator_refines_successful_resource_creation() {
    let input = r#"
program resource
  use mpi_f08, only: mpi_comm, mpi_comm_world, mpi_comm_dup, mpi_comm_free, mpi_abort, mpi_success
  implicit none (type, external)
  type(mpi_comm) :: comm
  integer :: ierr, ignored
  call mpi_comm_dup(mpi_comm_world, comm, ierror=ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ignored)
  call mpi_comm_free(comm, ierror=ierr)
  if (ierr /= mpi_success) call mpi_abort(mpi_comm_world, ierr, ignored)
end program
"#;
    let config = Config::default();
    let actual = diagnose(input, &config);
    assert!(
        !has(&actual, "OWN001") && !has(&actual, "PRL001"),
        "{actual:?}"
    );
    let actual = diagnose(
        &input.replace("call mpi_comm_free(comm, ierror=ierr)", "continue"),
        &config,
    );
    assert!(has(&actual, "OWN001"), "{actual:?}");
}

#[test]
fn renamed_calls_detect_uncompleted_requests_and_rank_dependent_collectives() {
    let input = r#"
program requests
  use mpi_f08, only: mpi_comm_world, mpi_integer, mpi_request, mpi_success, &
       find_rank => mpi_comm_rank, send => mpi_isend, mpi_barrier
  implicit none (type, external)
  type(mpi_request) :: request
  integer :: ierr, rank, item
  call find_rank(mpi_comm_world, rank, ierror=ierr)
  if (ierr /= mpi_success) error stop
  call send(item, 1, mpi_integer, 0, 12, mpi_comm_world, request, ierror=ierr)
  if (ierr /= mpi_success) error stop
  if (rank == 0) then
    call mpi_barrier(mpi_comm_world, ierror=ierr)
    if (ierr /= mpi_success) error stop
  end if
end program
"#;
    let actual = diagnose(input, &Config::default());
    assert!(
        has(&actual, "PRL003") && has(&actual, "PRL005"),
        "{actual:?}"
    );
    assert!(!has(&actual, "PRL001"), "{actual:?}");
}

#[test]
fn mpi_lookalike_project_procedure_is_not_a_library_call() {
    let input = r#"
program application
  implicit none (type, external)
  call mpi_barrier()
contains
  subroutine mpi_barrier()
    implicit none (type, external)
  end subroutine
end program
"#;
    assert!(diagnose(input, &Config::default()).is_empty());
}

#[test]
fn heap_assignment_and_temporary_detection_respect_sections() {
    let input = r#"
module arrays
  implicit none (type, external)
contains
  subroutine caller(input)
    real, intent(in) :: input(4, 4)
    real, allocatable :: values(:)
    values = [1.0, 2.0, 3.0, 4.0]
    values(:) = [1.0, 2.0, 3.0, 4.0]
    call consume(input(2, :))
    call consume(input(:, 2))
    call consume(input(::-1, 2))
    call consume(arg=input(4:1:-1, 2))
    call consume_view(input(1:4:2, 2))
  end subroutine
  subroutine consume(arg)
    real, contiguous, intent(in) :: arg(:)
  end subroutine
  subroutine consume_view(arg)
    real, intent(in) :: arg(:)
  end subroutine
end module
"#;
    let diagnostics = diagnose(input, &Config::default());
    assert_eq!(
        diagnostics
            .iter()
            .filter(|item| item.code == "MEM004")
            .count(),
        1,
        "{diagnostics:?}"
    );
    let mut packed: Vec<_> = diagnostics
        .iter()
        .filter(|item| item.code == "MEM005")
        .map(|item| item.location.line)
        .collect();
    packed.sort_unstable();
    assert_eq!(packed, [10, 13], "{diagnostics:?}");
}

#[test]
fn pointer_checks_distinguish_conditional_and_null_associations() {
    let input = r#"
subroutine check_pointer(condition)
  use, intrinsic :: iso_c_binding, only: bridge => c_f_pointer, c_loc
  implicit none (type, external)
  logical, intent(in) :: condition
  real, target :: values(3)
  real, pointer :: guarded(:), uncertain(:)
  nullify(guarded)
  if (associated(guarded)) print *, guarded(1)
  call bridge(c_loc(values), guarded, [3])
  if (associated(guarded)) print *, guarded(1)
  if (condition) uncertain => values
  print *, uncertain(1)
end subroutine
"#;
    let diagnostics = diagnose(input, &Config::default());
    let expected: Vec<_> = diagnostics
        .iter()
        .map(|item| (item.code.as_str(), item.severity))
        .collect();
    assert_eq!(expected, [("MEM002", Severity::Warning)], "{diagnostics:?}");
}

#[test]
fn petsc_views_are_restored_using_the_same_owner_and_mode() {
    let input = r#"
subroutine read_view(vector, other)
  use petscvec, only: tvec, vecgetarrayread, vecrestorearrayread, vecrestorearray
  implicit none (type, external)
  type(tvec), intent(in) :: vector, other
  real, pointer :: view(:)
  integer :: ierr
  call vecgetarrayread(vector, view, ierr)
  if (ierr /= 0) return
  print *, view(1)
  call vecrestorearrayread(vector, view, ierr)
  if (ierr /= 0) error stop
end subroutine
"#;
    let config = Config::default();
    assert!(diagnose(input, &config).is_empty());
    let mode = input.replace(
        "call vecrestorearrayread(vector",
        "call vecrestorearray(vector",
    );
    let wrong_mode = diagnose(&mode, &config);
    assert!(
        wrong_mode
            .iter()
            .any(|item| item.code == "OWN002" && item.severity == Severity::Error),
        "{wrong_mode:?}"
    );
    let owner = input.replace(
        "call vecrestorearrayread(vector",
        "call vecrestorearrayread(other",
    );
    let wrong_owner = diagnose(&owner, &config);
    assert!(
        wrong_owner
            .iter()
            .any(|item| item.code == "OWN002" && item.severity == Severity::Warning),
        "{wrong_owner:?}"
    );
}

#[test]
fn handles_in_types_require_cleanup_and_receive_copy_diagnostics() {
    let input = r#"
module holders
  use mpi_f08, only: mpi_comm, mpi_comm_dup, mpi_comm_free
  implicit none (type, external)
  private
  type :: holder
    private
    type(mpi_comm) :: owned
  contains
    final :: close_holder
  end type
contains
  subroutine open_holder(self, comm)
    type(holder), intent(inout) :: self
    type(mpi_comm), intent(in) :: comm
    call mpi_comm_dup(comm, self%owned)
  end subroutine
  subroutine close_holder(self)
    type(holder), intent(inout) :: self
    call mpi_comm_free(self%owned)
  end subroutine
  subroutine copy_holder(destination, source)
    type(holder), intent(inout) :: destination
    type(holder), intent(in) :: source
    destination = source
  end subroutine
end module
"#;
    let config = Config::default();
    let good_cleanup = diagnose(input, &config);
    assert!(
        has(&good_cleanup, "OWN003") && !has(&good_cleanup, "OWN004"),
        "{good_cleanup:?}"
    );
    let bad_cleanup = input.replace("call mpi_comm_free(self%owned)", "continue");
    let missing = diagnose(&bad_cleanup, &config);
    assert!(
        missing
            .iter()
            .any(|item| item.code == "OWN004" && item.severity == Severity::Warning),
        "{missing:?}"
    );
    let invalid = input.replace(
        "subroutine close_holder(self)\n    type(holder)",
        "subroutine close_holder(self)\n    class(holder)",
    );
    let invalid = diagnose(&invalid, &config);
    assert!(
        invalid
            .iter()
            .any(|item| item.code == "OWN004" && item.severity == Severity::Error),
        "{invalid:?}"
    );
    let repeated = input.replace(
        "call mpi_comm_dup(comm, self%owned)",
        "call mpi_comm_dup(comm, self%owned)\n    call mpi_comm_dup(comm, self%owned)",
    );
    let overwritten = diagnose(&repeated, &config);
    assert!(has(&overwritten, "OWN005"), "{overwritten:?}");
}

#[test]
fn mpi_status_carriers_are_tracked_across_control_flow_joins() {
    let input = r#"
subroutine trial(retry, result_code)
  use mpi_f08, only: mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_success
  implicit none (type, external)
  logical, intent(in) :: retry
  integer, intent(inout) :: result_code
  integer :: ierr, rank, saved, count
  saved = 0
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  saved = ierr
  if (saved /= mpi_success) return
  call mpi_comm_size(mpi_comm_world, count, ierror=result_code)
end subroutine
"#;
    let config = Config::default();
    assert!(diagnose(input, &config).is_empty());
    let conditional = input.replace("saved = ierr", "if (retry) saved = ierr");
    let examples = diagnose(&conditional, &config);
    assert!(
        examples
            .iter()
            .any(|item| item.code == "PRL001" && item.severity == Severity::Warning),
        "{examples:?}"
    );
    let omitted = input.replace(
        "mpi_comm_world, count, ierror=result_code",
        "mpi_comm_world, count",
    );
    let examples = diagnose(&omitted, &config);
    assert!(
        examples
            .iter()
            .any(|item| item.code == "PRL001" && item.severity == Severity::Error),
        "{examples:?}"
    );
}

#[test]
fn split_petsc_collective_tracks_status_and_argument_order() {
    let input = r#"
subroutine scatter()
  use petscvec, only: vecscatterbegin, vecscatterend, insert_values, scatter_forward
  implicit none (type, external)
  integer :: scatter_context, field, alternative, ierr
  call vecscatterbegin(scatter_context, field, alternative, insert_values, scatter_forward, ierr)
  if (ierr /= 0) return
  call vecscatterend(scatter_context, field, alternative, insert_values, scatter_forward, ierr)
  if (ierr /= 0) error stop
end subroutine
"#;
    let config = Config::default();
    let diagnostics = diagnose(input, &config);
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let bad = input.replace(
        "vecscatterend(scatter_context, field, alternative, insert_values, scatter_forward, ierr)",
        "vecscatterend(scatter_context, alternative, field, insert_values, scatter_forward, ierr)",
    );
    let diagnostics = diagnose(&bad, &config);
    assert!(has(&diagnostics, "PRL003"), "{diagnostics:?}");
    let dropped = input.replace("  if (ierr /= 0) error stop\n", "");
    let diagnostics = diagnose(&dropped, &config);
    assert!(has(&diagnostics, "PRL002"), "{diagnostics:?}");
}

#[test]
fn mpi_test_flag_proves_completion_on_the_success_path_only() {
    let input = r#"
subroutine communicate()
  use mpi_f08, only: mpi_integer, mpi_comm_world, mpi_barrier, mpi_ibARRIER, mpi_test, mpi_request_free, mpi_success
  implicit none (type, external)
  logical :: done, complete
  integer :: payload, request, ierr, status
  call MPI_Ibarrier(mpi_comm_world, request, ierr)
  if (ierr /= mpi_success) return
  call MPI_Test(request, done, status, ierr)
  if (ierr /= mpi_success) error stop
  complete = done
  if (.not. complete) then
    call mpi_test(request, done, status, ierr)
    if (ierr /= mpi_success) error stop
    if (.not. done) error stop
  end if
end subroutine
"#;
    let config = Config::default();
    let diagnostics = diagnose(input, &config);
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let pending = input.replace("    if (.not. done) error stop\n", "");
    let diagnostics = diagnose(&pending, &config);
    assert!(has(&diagnostics, "PRL003"), "{diagnostics:?}");
    let costly = input.replace(
        "  call MPI_Test(request, done, status, ierr)",
        "  call mpi_request_free(request, ierr)",
    );
    let diagnostics = diagnose(&costly, &config);
    assert!(
        diagnostics
            .iter()
            .any(|item| item.code == "PRL004" && item.severity == Severity::Error),
        "{diagnostics:?}"
    );
}

#[test]
fn explicit_parallel_regions_track_shared_static_writes_and_localities() {
    let input = r#"
module records
  implicit none (type, external)
  type :: holder
    integer :: component
  end type
  type(holder) :: shared_record
  integer :: indexed(2), counter, independent, racing
contains
  subroutine compute()
    integer :: index, temporary
    integer, save :: cached
    associate (alias => shared_record%component)
      do concurrent (index=1:2) local(temporary)
        indexed(index) = index
        alias = index
        cached = index
        temporary = index
      end do
    end associate
    !$omp parallel do private(independent)
    do index=1,2
      !$omp atomic update
      counter = counter + 1
      independent = index
      racing = index
    end do
    !$omp end parallel do
  end subroutine
end module
"#;
    let mut config = Config::default();
    let diagnostics = diagnose(input, &config);
    let effects: Vec<_> = diagnostics
        .iter()
        .filter(|item| item.code == "PRL006")
        .map(|item| item.severity)
        .collect();
    assert_eq!(
        effects,
        [Severity::Error, Severity::Error, Severity::Warning],
        "{diagnostics:?}"
    );
    config.safety.check_do_concurrent = false;
    config.safety.check_openmp = false;
    assert!(diagnose(input, &config).is_empty());
}

#[test]
fn separate_pure_procedure_writes_its_result_and_dummy_locally() {
    let input = r#"
module owner
    implicit none (type, external)
    integer :: shared
    interface
        module pure function compute(input) result(value)
            integer, intent(in) :: input
            integer :: value
        end function
        module pure subroutine initialize(value)
            integer, intent(out) :: value
        end subroutine
    end interface
end module
submodule (owner) implementation
contains
    module procedure compute
        value = input
        shared = input
    end procedure
    module procedure initialize
        value = 0
    end procedure
end submodule
"#;
    let diagnostics = diagnose(input, &Config::default());
    let effects: Vec<_> = diagnostics
        .iter()
        .filter(|item| item.code == "PRL006")
        .collect();
    assert_eq!(effects.len(), 1, "{diagnostics:?}");
    assert_eq!(effects[0].severity, Severity::Error);
    assert!(effects[0].message.contains("'shared'"));
}
