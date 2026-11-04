program gibbs_derivatives
    use iso_fortran_env, only: real64, int64
    use ieee_arithmetic, only: ieee_is_finite
    use iapws, only: type_iapws97
    implicit none
    type(type_iapws97) :: water
    real(real64), parameter :: temperatures(9) = [300.d0, 500.d0, 600.d0, &
        400.d0, 700.d0, 1000.d0, 1100.d0, 1500.d0, 2200.d0]
    real(real64), parameter :: pressures(9) = [3.d6, 10.d6, 30.d6, &
        1.d5, 3.d6, 20.d6, 1.d6, 30.d6, 40.d6]
    real(real64) :: dt(9), dp(9), cp_plus(9), cp_minus(9), ref(9), ht, hp(9)
    real(real64) :: bulk_dt(9), bulk_dp(9)
    integer :: i, step
    call water%initialize()
    call water%calc_dcp_dT(temperatures, pressures, dt)
    call water%calc_dcp_dp(temperatures, pressures, dp)
    call water%calc_cp_derivatives(temperatures,pressures,bulk_dt,bulk_dp)
    if(any(transfer(dt,[0_int64])/=transfer(bulk_dt,[0_int64])).or. &
       any(transfer(dp,[0_int64])/=transfer(bulk_dp,[0_int64]))) error stop 'bulk cp derivative compatibility'
    dt=bulk_dt; dp=bulk_dp
    ! Two independent steps catch cancellation and the Pa versus MPa conversion.
    do step = 1, 2
        ht = 0.01d0 / real(step, real64)
        hp = pressures * (1.d-4 / real(step, real64))
        call water%calc_cp(temperatures + ht, pressures, cp_plus)
        call water%calc_cp(temperatures - ht, pressures, cp_minus)
        ref = (cp_plus - cp_minus) / (2.d0 * ht)
        do i = 1, size(dt)
            call check(dt(i), ref(i), 1.d-7, 'temperature', i)
        end do
        call water%calc_cp(temperatures, pressures + hp, cp_plus)
        call water%calc_cp(temperatures, pressures - hp, cp_minus)
        ref = (cp_plus - cp_minus) / (2.d0 * hp)
        do i = 1, size(dp)
            call check(dp(i), ref(i), 1.d-12, 'pressure', i)
        end do
    end do
contains
    subroutine check(actual, expected, atol, label, index)
        real(real64), intent(in) :: actual, expected, atol
        character(*), intent(in) :: label
        integer, intent(in) :: index
        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) error stop 'nonfinite derivative'
        if (abs(actual - expected) > atol + 1.d-5 * abs(expected)) then
            print *, label, index, actual, expected
            error stop 'Gibbs specific heat derivative mismatch'
        end if
    end subroutine check
end program gibbs_derivatives
