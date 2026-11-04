program ice_derivatives
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic
    use :: module_iapws06_iceIh, only: type_iapws06
    implicit none
    type(type_iapws06) :: ice
    real(real64), parameter :: temperatures(5) = [20d0, 100d0, 200d0, 250d0, 273d0]
    real(real64), parameter :: pressures(5) = [0d0, 1d5, 5d7, 1d8, 2.1d8]
    real(real64) :: dt(5), dp(5), plus, minus, expected, step, bad(8), values(8), pbad(8)
    integer :: i, j
    call ice%initialize()
    call ice%calc_dcp_dT(temperatures, pressures, dt)
    call ice%calc_dcp_dP(temperatures, pressures, dp)
    do i=1,5
        do j=1,3
            step = 0.01d0 / real(j,real64)
            call ice%calc_cp(temperatures(i)+step, pressures(i), plus)
            call ice%calc_cp(temperatures(i)-step, pressures(i), minus)
            expected=(plus-minus)/(2d0*step)
            call check(dt(i), expected, 1d-5)
            step=1000d0/real(j,real64)
            call ice%calc_cp(temperatures(i), pressures(i)+step, plus)
            call ice%calc_cp(temperatures(i), pressures(i)-step, minus)
            expected=(plus-minus)/(2d0*step)
            call check(dp(i), expected, 1d-5)
        end do
    end do
    bad=[0d0,-1d0,273.17d0,200d0,200d0,ieee_value(0d0,ieee_quiet_nan), &
         ieee_value(0d0,ieee_positive_inf),200d0]
    pbad=[0d0,0d0,0d0,-1d0,210000001d0,0d0,0d0,ieee_value(0d0,ieee_positive_inf)]
    call ice%calc_dcp_dT(bad,pbad,values)
    if (.not.all(ieee_is_nan(values))) error stop 'invalid ice dcp/dT must be NaN'
    call ice%calc_dcp_dP(bad,pbad,values)
    if (.not.all(ieee_is_nan(values))) error stop 'invalid ice dcp/dp must be NaN'
    call ice%calc_dcp_dT(273.16d0,0d0,plus)
    call ice%calc_dcp_dP(273.16d0,210000000d0,minus)
    if (.not.ieee_is_finite(plus) .or. .not.ieee_is_finite(minus)) error stop 'valid domain boundaries'
    print *, 'Ice analytic derivatives passed'
contains
    subroutine check(actual,reference,tolerance)
        real(real64), intent(in) :: actual,reference,tolerance
        if (.not.ieee_is_finite(actual)) error stop 'nonfinite ice derivative'
        if (abs(actual-reference)>tolerance*max(abs(reference),1d-12)) then
            print *, actual,reference
            error stop 'ice derivative mismatch'
        end if
    end subroutine
end program
