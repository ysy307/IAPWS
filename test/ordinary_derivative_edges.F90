program ordinary_derivative_edges
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic
    use iapws, only: type_iapws97
    implicit none
    type(type_iapws97) :: water
    real(real64) :: ts(20), ps(20), dt(20), dp(20), nan, inf
    real(real64), parameter :: valid_t(4) = [300d0, 700d0, 650d0, 1500d0]
    real(real64), parameter :: valid_p(4) = [3d6, 3d6, 25d6, 30d6]
    real(real64) :: boundary_t(3), boundary_p(3), dt_bound(3), dp_bound(3)
    integer :: i

    call water%initialize()
    nan = ieee_value(0d0, ieee_quiet_nan)
    inf = ieee_value(0d0, ieee_positive_inf)
    ts = [valid_t, 0d0, -1d0, 243.15d0-1d-6, 2273.15d0+1d-6, &
          300d0, 300d0, 300d0, 1500d0, nan, inf, -inf, 300d0, &
          300d0, 300d0, 647.096d0, 647.096d0]
    ps = [valid_p, 1d5, 1d5, 1d5, 1d5, 0d0, -1d0, 100d6+1d0, &
          50d6+1d0, 1d5, 1d5, 1d5, nan, inf, -inf, 22.064d6, 22.064d6]
    ! A pure procedure and mixed elemental arrays must retain valid elements
    ! while invalid elements and the critical singularity independently yield NaN.
    call evaluate_pure(water, ts, ps, dt, dp)
    if (.not.all(ieee_is_finite(dt(1:4)))) error stop 'valid mixed-array dcp/dT'
    if (.not.all(ieee_is_finite(dp(1:4)))) error stop 'valid mixed-array dcp/dp'
    do i = 5, size(ts)
        if (.not.ieee_is_nan(dt(i)) .or. .not.ieee_is_nan(dp(i))) then
            print *, 'invalid state', i, ts(i), ps(i), dt(i), dp(i)
            error stop 'invalid ordinary derivative must yield NaN'
        end if
    end do

    ! Preserve the existing ordinary-property lower bound and include exact
    ! region-2 and region-5 upper bounds (pressure derivatives remain per Pa).
    boundary_t = [243.15d0, 1073.15d0, 2273.15d0]
    boundary_p = [1d5, 100d6, 50d6]
    call evaluate_pure(water, boundary_t, boundary_p, dt_bound, dp_bound)
    if (.not.all(ieee_is_finite(dt_bound))) error stop 'valid boundary dcp/dT'
    if (.not.all(ieee_is_finite(dp_bound))) error stop 'valid boundary dcp/dp'

    ! At a fit join, differences must stay on the branch selected at the join.
    ! Region 1 owns T=623.15 K at 20 MPa; region 2 owns T=1073.15 K.
    call check_join(623.15d0, 20d6, -1d0)
    call check_join(1073.15d0, 30d6, -1d0)
    call check_join(623.15d0+0.01d0, 20d6, 1d0)
    call check_join(1073.15d0+0.01d0, 30d6, 1d0)
    print *, 'Ordinary derivative invalid inputs, pure arrays, bounds and joins passed'
contains
    pure subroutine evaluate_pure(model, temperature, pressure, dT, dP)
        type(type_iapws97), intent(in) :: model
        real(real64), intent(in) :: temperature(:), pressure(:)
        real(real64), intent(out) :: dT(:), dP(:)
        call model%calc_dcp_dT(temperature, pressure, dT)
        call model%calc_dcp_dp(temperature, pressure, dP)
    end subroutine

    subroutine check_join(temperature, pressure, direction)
        real(real64), intent(in) :: temperature, pressure, direction
        real(real64) :: actual, expected, c0, c1, c2, h, cp_plus, cp_minus
        integer :: step
        call water%calc_dcp_dT(temperature, pressure, actual)
        do step = 1, 2
            h = direction * 1d-3 / real(step, real64)
            call water%calc_cp(temperature, pressure, c0)
            call water%calc_cp(temperature+h, pressure, c1)
            call water%calc_cp(temperature+2*h, pressure, c2)
            expected = (-3*c0+4*c1-c2)/(2*h)
            call close(actual, expected, 'branch-local temperature')
        end do
        call water%calc_dcp_dp(temperature, pressure, actual)
        do step = 1, 2
            h = 100d0 / real(step, real64)
            call water%calc_cp(temperature, pressure+h, cp_plus)
            call water%calc_cp(temperature, pressure-h, cp_minus)
            expected = (cp_plus-cp_minus)/(2*h)
            call close(actual, expected, 'join pressure')
        end do
    end subroutine

    subroutine close(actual, expected, label)
        real(real64), intent(in) :: actual, expected
        character(*), intent(in) :: label
        if (.not.ieee_is_finite(actual) .or. .not.ieee_is_finite(expected)) error stop 'nonfinite join derivative'
        if (abs(actual-expected)>1d-4*max(abs(expected),1d-10)) then
            print *, label, actual, expected
            error stop 'join derivative mismatch'
        end if
    end subroutine
end program
