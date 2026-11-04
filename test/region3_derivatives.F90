program region3_derivatives
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iapws, only: type_iapws97
    implicit none
    type(type_iapws97) :: water
    real(real64), parameter :: temperatures(9) = [650d0,700d0,640d0,630d0,627.15d0,631.15d0,635.15d0,639.15d0,646.15d0]
    real(real64), parameter :: pressures(9) = [25d6,30d6,20d6,20d6,17d6,18d6,19d6,20d6,22d6]
    real(real64) :: dt(9), dp(9), plus, minus, expected, ht, hp
    integer :: i, step
    call water%initialize()
    ! Elemental array calls also exercise the public dispatch and density solve.
    call water%calc_dcp_dT(temperatures, pressures, dt)
    call water%calc_dcp_dp(temperatures, pressures, dp)
    do i = 1, size(temperatures)
        do step = 1, 2
            ht = 0.002_real64 / real(step, real64)
            hp = 100.0_real64 / real(step, real64)
            call water%calc_cp(temperatures(i)+ht, pressures(i), plus)
            call water%calc_cp(temperatures(i)-ht, pressures(i), minus)
            expected = (plus-minus)/(2*ht)
            call check(dt(i), expected, 'temperature derivative', i)
            call water%calc_cp(temperatures(i), pressures(i)+hp, plus)
            call water%calc_cp(temperatures(i), pressures(i)-hp, minus)
            expected = (plus-minus)/(2*hp)
            call check(dp(i), expected, 'pressure derivative', i)
        end do
    end do
contains
    subroutine check(actual, expected, label, state)
        real(real64), intent(in) :: actual, expected
        character(*), intent(in) :: label
        integer, intent(in) :: state
        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) error stop 'nonfinite derivative'
        if (abs(actual-expected) > 1.e-4_real64*max(abs(expected),1.e-9_real64)) then
            print *, label, state, actual, expected
            error stop 'region 3 derivative mismatch'
        end if
    end subroutine check
end program region3_derivatives
