program main
    use, intrinsic :: iso_fortran_env
    use :: module_iapws
    implicit none

    call test_iapws95()
    call test_iapws97()

contains
    subroutine test_iapws95()
        implicit none
        integer(int32) :: unit
        integer(int32), parameter :: test_case = 11
        real(real64) :: T(test_case), rho(test_case)
        real(real64) :: p(test_case), cv(test_case), w(test_case), s(test_case)
        integer(int32) :: i
        logical :: exists

        type(type_iapws95) :: iapws95_model

        type(type_iapws_property) :: props(test_case), prop

        open (unit=10, file="/workspaces/IAPWS/validations/iapws95/Helmholtz.dat", status="old", action="read")
        read (10, *)
        do i = 1, test_case
            read (10, *) T(i), rho(i), p(i), cv(i), w(i), s(i)
        end do
        close (10)

        inquire (file="/workspaces/IAPWS/validations/iapws95/test_iapws95.log", exist=exists)
        if (.not. exists) then
            open (newunit=unit, file="/workspaces/IAPWS/validations/iapws95/test_iapws95.log", &
                  status="new", action="write")
        else
            open (newunit=unit, file="/workspaces/IAPWS/validations/iapws95/test_iapws95.log", &
                  status="old", action="write", position="append")
        end if

        write (unit, '(a)') "# IAPWS-95 Validation Test Log ("//trim(get_compiler_name())//")"

        ! Calculate ideal Helmholtz properties
        call iapws95_model%initialize()
        call iapws95_model%calc_properties(T(:), rho(:), props(:))

        call check_variables(unit, [(props(i)%p * 1.0d-6, i=1, test_case)], p, T(:), "Pressure", [(i, i=1, test_case)])
        call check_variables(unit, [(props(i)%cv * 1.0d-3, i=1, test_case)], cv, T(:), "cv", [(i, i=1, test_case)])
        call check_variables(unit, [(props(i)%w, i=1, test_case)], w, T(:), "w", [(i, i=1, test_case)])
        call check_variables(unit, [(props(i)%s * 1.0d-3, i=1, test_case)], s, T(:), "s", [(i, i=1, test_case)])

        close (unit)

    end subroutine test_iapws95

    subroutine test_iapws97()
        implicit none
        integer(int32) :: unit
        integer(int32), parameter :: test_gibbs = 9
        real(real64) :: T_g(test_gibbs), p_g(test_gibbs)
        real(real64) :: nu_g_exa(test_gibbs)
        real(real64) :: u_g_exa(test_gibbs), h_g_exa(test_gibbs)
        real(real64) :: s_g_exa(test_gibbs), cp_g_exa(test_gibbs)
        real(real64) :: w_g_exa(test_gibbs)
        integer(int32) :: region_gibbs(test_gibbs)
        integer(int32), parameter :: test_helmholtz = 3
        real(real64) :: T_h(test_helmholtz), rho_h(test_helmholtz)
        real(real64) :: p_h_exa(test_helmholtz), u_h_exa(test_helmholtz)
        real(real64) :: h_h_exa(test_helmholtz), s_h_exa(test_helmholtz)
        real(real64) :: cp_h_exa(test_helmholtz), w_h_exa(test_helmholtz)
        integer(int32) :: region_h(test_helmholtz)
        integer(int32) :: i
        logical :: exists

        type(type_iapws97) :: iapws97_model
        type(type_iapws_property) :: props_g(test_gibbs), props_h(test_helmholtz)
        type(type_iapws_gamma_property) :: prop_g
        type(type_iapws_phi_property) :: prop_h

        open (unit=10, file="/workspaces/IAPWS/validations/iapws97/test_iapws97.dat", status="old", action="read")
        read (10, *)
        do i = 1, test_gibbs
            read (10, *) region_gibbs(i), T_g(i), p_g(i), nu_g_exa(i), h_g_exa(i), u_g_exa(i), &
                s_g_exa(i), cp_g_exa(i), w_g_exa(i)
        end do
        read (10, *)
        read (10, *)
        do i = 1, test_helmholtz
            read (10, *) region_h(i), T_h(i), rho_h(i), p_h_exa(i), h_h_exa(i), u_h_exa(i), &
                s_h_exa(i), cp_h_exa(i), w_h_exa(i)
        end do
        close (10)

        inquire (file="/workspaces/IAPWS/validations/iapws97/test_iapws97.log", exist=exists)
        if (.not. exists) then
            open (newunit=unit, file="/workspaces/IAPWS/validations/iapws97/test_iapws97.log", &
                  status="new", action="write")
        else
            open (newunit=unit, file="/workspaces/IAPWS/validations/iapws97/test_iapws97.log", &
                  status="old", action="write", position="append")
        end if

        write (unit, '(a)') "# IAPWS-97 Validation Test Log ("//trim(get_compiler_name())//")"

        call iapws97_model%initialize()
        call iapws97_model%region1%calc_properties(T_g(1:3), p_g(1:3), props_g(1:3))
        call iapws97_model%region2%calc_properties(T_g(4:6), p_g(4:6), props_g(4:6))
        call iapws97_model%region3%calc_properties(T_h(:), rho_h(:), props_h(:))
        call iapws97_model%region5%calc_properties(T_g(7:9), p_g(7:9), props_g(7:9))

        write (unit, '(a)') "## Gibbs Formulation Tests"
        call check_variables(unit, [(props_g(i)%nu, i=1, test_gibbs)], &
                             nu_g_exa, T_g(:), "Specific Volume", [(i, i=1, test_gibbs)])
        call check_variables(unit, [(props_g(i)%u, i=1, test_gibbs)], &
                             u_g_exa, T_g(:), "Internal Energy", [(i, i=1, test_gibbs)])
        call check_variables(unit, [(props_g(i)%h, i=1, test_gibbs)], &
                             h_g_exa, T_g(:), "Enthalpy", [(i, i=1, test_gibbs)])
        call check_variables(unit, [(props_g(i)%s, i=1, test_gibbs)], &
                             s_g_exa, T_g(:), "Entropy", [(i, i=1, test_gibbs)])
        call check_variables(unit, [(props_g(i)%cp, i=1, test_gibbs)], &
                             cp_g_exa, T_g(:), "cp", [(i, i=1, test_gibbs)])
        call check_variables(unit, [(props_g(i)%w, i=1, test_gibbs)], &
                             w_g_exa, T_g(:), "w", [(i, i=1, test_gibbs)])

        write (unit, '(a)') "## Helmholtz Formulation Tests"
        call check_variables(unit, [(props_h(i)%p, i=1, test_helmholtz)], &
                             p_h_exa, T_h(:), "Pressure", [(i, i=1, test_helmholtz)])
        call check_variables(unit, [(props_h(i)%u, i=1, test_helmholtz)], &
                             u_h_exa, T_h(:), "Internal Energy", [(i, i=1, test_helmholtz)])
        call check_variables(unit, [(props_h(i)%h, i=1, test_helmholtz)], &
                             h_h_exa, T_h(:), "Enthalpy", [(i, i=1, test_helmholtz)])
        call check_variables(unit, [(props_h(i)%s, i=1, test_helmholtz)], &
                             s_h_exa, T_h(:), "Entropy", [(i, i=1, test_helmholtz)])
        call check_variables(unit, [(props_h(i)%cp, i=1, test_helmholtz)], &
                             cp_h_exa, T_h(:), "cp", [(i, i=1, test_helmholtz)])
        call check_variables(unit, [(props_h(i)%w, i=1, test_helmholtz)], &
                             w_h_exa, T_h(:), "w", [(i, i=1, test_helmholtz)])

    end subroutine test_iapws97

    function get_compiler_name() result(name)
        implicit none
        character(len=100) :: name

#ifdef __GFORTRAN__
        name = "GCC Fortran Compiler"
#elif defined(__INTEL_COMPILER)
        name = "Intel Fortran Compiler"
#elif defined(__PGI) || defined(__NVCOMPILER)
        name = "NVIDIA Fortran Compiler"
#else
        name = "Unknown Compiler"
#endif
    end function get_compiler_name

    subroutine check_variable(unit, v, v_exa, T_in, v_name, id)
        implicit none
        integer(int32), intent(in) :: unit
        real(real64), intent(in) :: v, v_exa
        real(real64), intent(in) :: T_in
        character(len=*), intent(in) :: v_name
        integer(int32), intent(in), optional :: id

        real(real64), parameter :: tol = 5.0d-9
        real(real64) :: rel_diff
        real(real64) :: current_tol

        if (abs(T_in - 647.0d0) < 1.0d0) then
            current_tol = 5.0d-6
        else
            current_tol = 1.0d-9 ! それ以外は厳しく
        end if

        if (abs(v_exa) > 0.0d0) then
            rel_diff = abs(v - v_exa) / abs(v_exa)
        else
            rel_diff = abs(v - v_exa)
        end if

        if (rel_diff > tol) then
            write (unit, '(a)') "**FAIL**: `"//v_name//"`"
            write (unit, '(a)') ""
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') "ID", "computed", "expected", "rel_diff"
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') &
                repeat('-', 6), repeat('-', 20), repeat('-', 20), repeat('-', 20)
            if (present(id)) then
                write (unit, '("|",i6,"|",es20.8,"|",es20.8,"|",es20.8,"|")') id, v, v_exa, rel_diff
            else
                write (unit, '("|",a6,"|",es20.8,"|",es20.8,"|",es20.8,"|")') "-", v, v_exa, rel_diff
            end if
            write (unit, '(a)') ""
        else
            write (unit, '(a)') "PASS: `"//v_name//"`"
            write (unit, '(a)') ""
        end if
    end subroutine check_variable

    subroutine check_variables(unit, v, v_exa, T_in, v_name, ids)
        implicit none
        integer(int32), intent(in) :: unit
        real(real64), intent(in) :: v(:)
        real(real64), intent(in) :: v_exa(:)
        real(real64), intent(in) :: T_in(:)
        character(len=*), intent(in) :: v_name
        integer(int32), intent(in), optional :: ids(:)

        real(real64), parameter :: tol = 5.0d-9

        real(real64), allocatable :: rel_diff(:)
        integer(int32) :: i, n
        integer(int32) :: fail_count
        real(real64), allocatable :: current_tol(:)

        n = size(v)
        allocate (rel_diff(n))
        allocate (current_tol(n))

        do i = 1, n
            if (abs(T_in(i) - 647.0d0) < 1.0d0) then
                current_tol(i) = 5.0d-6
            else
                current_tol(i) = 1.0d-9
            end if
        end do

        where (abs(v_exa) > 0.0d0)
            rel_diff = abs(v - v_exa) / abs(v_exa)
        elsewhere
            rel_diff = abs(v - v_exa)
        end where

        if (any(rel_diff > tol)) then
            fail_count = count(rel_diff > tol)
            write (unit, '(a,i0,a)') "**FAIL**: `"//v_name//"` (Failed count: ", fail_count, ")"
            write (unit, '(a)') ""
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') "ID", "computed", "expected", "rel_diff"
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') &
                repeat('-', 6), repeat('-', 20), repeat('-', 20), repeat('-', 20)

            do i = 1, n
                if (rel_diff(i) > tol) then
                    if (present(ids)) then
                        write (unit, '("|",i6,"|",es20.8,"|",es20.8,"|",es20.8,"|")') ids(i), v(i), v_exa(i), rel_diff(i)
                    else
                        write (unit, '("|",i6,"|",es20.8,"|",es20.8,"|",es20.8,"|")') i, v(i), v_exa(i), rel_diff(i)
                    end if
                end if
            end do
            write (unit, '(a)') ""
        else
            write (unit, '(a)') "PASS: `"//v_name//"`"
            write (unit, '(a)') ""
        end if
    end subroutine check_variables

end program main
