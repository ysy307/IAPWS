program main
    use, intrinsic :: iso_fortran_env
    use :: module_iapws
    implicit none

    call test_iapws95()

contains
    subroutine test_iapws95()
        implicit none
        integer(int32) :: unit
        integer(int32), parameter :: test_case = 6
        real(real64) :: T(test_case), rho(test_case)
        real(real64) :: phi0_exa(test_case), phi0_d_exa(test_case)
        real(real64) :: phi0_dd_exa(test_case), phi0_t_exa(test_case)
        real(real64) :: phi0_tt_exa(test_case), phi0_dt_exa(test_case)
        real(real64) :: phir_exa(test_case), phir_d_exa(test_case)
        real(real64) :: phir_dd_exa(test_case), phir_t_exa(test_case)
        real(real64) :: phir_tt_exa(test_case), phir_dt_exa(test_case)
        integer(int32) :: i
        logical :: exists

        type(type_iapws95) :: iapws95_model

        type(type_iapws_phi_property) :: phi_props(test_case)

        open (unit=10, file="/workspaces/IAPWS/validations/iapws95/Helmholtz.dat", status="old", action="read")
        read (10, *) ! ヘッダー行をスキップ
        do i = 1, test_case
            read (10, *) T(i), rho(i), &
                phi0_exa(i), phi0_d_exa(i), phi0_dd_exa(i), &
                phi0_t_exa(i), phi0_tt_exa(i), phi0_dt_exa(i), &
                phir_exa(i), phir_d_exa(i), phir_dd_exa(i), &
                phir_t_exa(i), phir_tt_exa(i), phir_dt_exa(i)
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

#ifdef __GFORTRAN__
        write (unit, '(a)') "# IAPWS-95 Validation Test Log (GCC Fortran Compiler)"
#elif defined(__INTEL_COMPILER)
        write (unit, '(a)') "# IAPWS-95 Validation Test Log (Intel Fortran Compiler)"
#elif defined(__PGI) || defined(__NVCOMPILER)
        write (unit, '(a)') "# IAPWS-95 Validation Test Log (NVIDIA Fortran Compiler)"
#else
        write (unit, '(a)') "# IAPWS-95 Validation Test Log (Unknown Compiler)"
#endif

        ! Calculate ideal Helmholtz properties
        call iapws95_model%initialize()
        call iapws95_model%calc_phi0(iapws95_model%T_c / T(:), rho(:) / iapws95_model%rho_c, phi_props(:))
        call iapws95_model%calc_phir(iapws95_model%T_c / T(:), rho(:) / iapws95_model%rho_c, phi_props(:))

        call check_variables(unit, [(phi_props(i)%phi0, i=1, test_case)], phi0_exa, "phi0", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phi0_d, i=1, test_case)], phi0_d_exa, "phi0_d", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phi0_dd, i=1, test_case)], phi0_dd_exa, "phi0_dd", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phi0_t, i=1, test_case)], phi0_t_exa, "phi0_t", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phi0_tt, i=1, test_case)], phi0_tt_exa, "phi0_tt", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phi0_dt, i=1, test_case)], phi0_dt_exa, "phi0_dt", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir, i=1, test_case)], phir_exa, "phir", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir_d, i=1, test_case)], phir_d_exa, "phir_d", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir_dd, i=1, test_case)], phir_dd_exa, "phir_dd", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir_t, i=1, test_case)], phir_t_exa, "phir_t", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir_tt, i=1, test_case)], phir_tt_exa, "phir_tt", [(i, i=1, test_case)])
        call check_variables(unit, [(phi_props(i)%phir_dt, i=1, test_case)], phir_dt_exa, "phir_dt", [(i, i=1, test_case)])

        close (unit)

    end subroutine test_iapws95

    subroutine check_variable(unit, v, v_exa, v_name, id)
        implicit none
        integer(int32), intent(in) :: unit
        real(real64), intent(in) :: v, v_exa
        character(len=*), intent(in) :: v_name
        integer(int32), intent(in), optional :: id
        ! --- 修正: 有効数字9桁を保証するために 1.0d-9 に設定 ---
        real(real64), parameter :: tol = 5.0d-9

        real(real64) :: rel_diff

        if (abs(v_exa) > 0.0d0) then
            rel_diff = abs(v - v_exa) / abs(v_exa)
        else
            ! 真値が0の場合は絶対誤差で評価
            rel_diff = abs(v - v_exa)
        end if

        if (rel_diff > tol) then
            write (unit, '(a)') "**FAIL**: `"//v_name//"`"
            write (unit, '(a)') ""
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') "ID", "computed", "expected", "rel_diff"
            write (unit, '("|",a6,"|",a20,"|",a20,"|",a20,"|")') &
                repeat('-', 6), repeat('-', 20), repeat('-', 20), repeat('-', 20)
            if (present(id)) then
                ! es20.810 は有効数字約11桁表示されるため、9桁の確認には十分です
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

    subroutine check_variables(unit, v, v_exa, v_name, ids)
        implicit none
        integer(int32), intent(in) :: unit
        real(real64), intent(in) :: v(:)
        real(real64), intent(in) :: v_exa(:)
        character(len=*), intent(in) :: v_name
        integer(int32), intent(in), optional :: ids(:)

        ! --- 修正: 有効数字9桁用 ---
        real(real64), parameter :: tol = 5.0d-9

        real(real64), allocatable :: rel_diff(:)
        integer(int32) :: i, n
        integer(int32) :: fail_count

        n = size(v)
        allocate (rel_diff(n))

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
                ! --- 修正: 失敗した要素だけを出力してログを見やすくする ---
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
