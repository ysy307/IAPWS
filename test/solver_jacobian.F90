program solver_jacobian
    use iapws
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    type(type_water) :: water
    type(type_water_state) :: q, a, b, recovered
    type(type_water_jacobian) :: jac
    type(type_water_saturation) :: sat
    type(type_water_state) :: cells(3), parallel_cells(96)
    real(real64) :: t, p, dh, h, s, nan
    integer :: i

    do i = 1, 3
        water = new_water(merge(IF97, IAPWS95, i == 1))
        t = merge(260d0, 400d0, i == 3)
        p = 3d6
        call water%state_tp(t, p, q)
        call q%jacobian(jac)
        if (jac%status /= WATER_OK) error stop 'Jacobian status'
        if (abs(jac%ph(JAC_T,2) * q%properties%cp - 1d0) > 8d-15) error stop 'dT/dh identity'
        dh = 0.5d0
        call water%state_ph(p, q%properties%h + dh, a, guess_T=t)
        call water%state_ph(p, q%properties%h - dh, b, guess_T=t)
        if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'warm inverse'
        if (abs((a%properties%rho-b%properties%rho)/(2d0*dh)-jac%ph(JAC_RHO,2)) > 2d-8) &
            error stop 'density enthalpy derivative'
        call water%state_ps(p, q%properties%s, recovered, guess_T=t)
        if (recovered%status /= WATER_OK .or. abs(recovered%properties%T-t) > 2d-9) error stop 'warm entropy inverse'
        if (recovered%iterations /= 1) error stop 'exact warm state needs one evaluation'
        call water%state_tp(t, p, recovered, guess_rho=q%properties%rho)
        if (recovered%status /= WATER_OK .or. abs(recovered%properties%rho/q%properties%rho-1d0)>1d-12) &
            error stop 'warm density branch'
        call water%state_ph_checked(-1d0, q%properties%h, recovered, guess_T=t)
        if (recovered%status == WATER_OK .or. .not. ieee_is_finite(recovered%properties%T)) error stop 'invalid p,h'
    end do
    call water%state_tp_checked(-1d0, -1d0, q)
    if (q%status == WATER_OK .or. .not. ieee_is_finite(q%properties%rho)) error stop 'finite failure'
    do i = 1, 2
        water = new_water(merge(IF97, IAPWS95, i == 1))
        p = 1d6
        call water%saturation_p(p, sat)
        if (sat%status /= WATER_OK) error stop 'saturation status'
        h = 0.5d0 * (sat%vapor%h + sat%liquid%h)
        s = 0.5d0 * (sat%vapor%s + sat%liquid%s)
        call verify_mixture(LIQUID_VAPOR, 100d0)
    end do

    p = melting_pressure(270d0)
    call water%state_tp(270d0, p, a, phase=ICE_IH)
    call water%state_tp(270d0, p, b, phase=LIQUID)
    if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'melting endpoints'
    h = 0.5d0 * (a%properties%h + b%properties%h)
    s = 0.5d0 * (a%properties%s + b%properties%s)
    call verify_mixture(ICE_LIQUID, 1000d0)

    p = sublimation_pressure(240d0)
    call water%state_tp(240d0, p, a, phase=ICE_IH)
    call water%state_tp(240d0, p, b, phase=VAPOR)
    if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'sublimation endpoints'
    h = 0.5d0 * (a%properties%h + b%properties%h)
    s = 0.5d0 * (a%properties%s + b%properties%s)
    call verify_mixture(ICE_VAPOR, 0.001d0)

    call water%state_tp_checked([400d0, -1d0, 260d0], 3d6, cells)
    if (any(cells%status /= [WATER_OK, WATER_INVALID, WATER_OK])) error stop 'elemental status'
    if (.not. all(ieee_is_finite(cells%properties%rho))) error stop 'elemental finite values'

    nan = ieee_value(0d0, ieee_quiet_nan)
    call water%state_tp_checked(nan, 3d6, q)
    if (q%status /= WATER_INVALID .or. .not. ieee_is_finite(q%properties%T)) error stop 'nonfinite input'
    call water%state_tp_checked(235.15d0, 101325d0, q, phase=LIQUID, allow_metastable=.true.)
    call q%jacobian(jac)
    if (q%formulation /= G12_15 .or. jac%status /= WATER_METASTABLE .or. .not. jac%tp_available) &
        error stop 'supercooled Jacobian'
    water = new_water(IF97)
    call water%state_tp_checked(450d0, 1d6, q, phase=VAPOR, allow_metastable=.true.)
    call q%jacobian(jac)
    if (q%formulation /= IF97_METASTABLE_VAPOR .or. jac%status /= WATER_METASTABLE) &
        error stop 'metastable steam Jacobian'
    call water%state_tp_checked(650d0, 25d6, q)
    if (q%status /= WATER_OK .or. q%properties%region_id /= 3) error stop 'region 3 state'
    call water%state_tp_checked(650d0, 25d6, a, guess_rho=q%properties%rho)
    call water%state_tp_checked(650d0, 25d6, b, guess_rho=nan)
    if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'region 3 density guesses'
    if (max(abs(a%properties%rho-q%properties%rho), abs(b%properties%rho-q%properties%rho)) > 1d-8) &
        error stop 'region 3 warm density mismatch'
    call water%state_ph_checked(25d6, q%properties%h, a, region=3, guess_T=650d0)
    if (a%status /= WATER_OK .or. a%iterations /= 1) error stop 'region 3 warm enthalpy'
    !$omp parallel do default(none) shared(water,q,parallel_cells)
    do i = 1, size(parallel_cells)
        call water%state_ph_checked(merge(25d6, -1d0, mod(i,2) == 0), q%properties%h, parallel_cells(i), &
                                    region=3, guess_T=650d0)
    end do
    !$omp end parallel do
    if (any(parallel_cells(1::2)%status /= WATER_INVALID) .or. &
        any(parallel_cells(2::2)%status /= WATER_OK)) error stop 'parallel status isolation'
    if (any(parallel_cells(2::2)%properties%rho /= a%properties%rho)) error stop 'parallel result mismatch'

contains

    subroutine verify_mixture(phase, dp)
        implicit none
        integer, intent(in) :: phase
        real(real64), intent(in) :: dp
        call water%state_ph_checked(p, h, q, jac=jac)
        if (q%status /= WATER_OK .or. q%phase /= phase .or. jac%status /= WATER_OK) error stop 'mixture'
        if (jac%tp_available .or. .not. jac%ph_available .or. .not. jac%ps_available) error stop 'mixture coordinates'
        if (q%properties%cp /= 0d0 .or. jac%ph(JAC_T,2) /= 0d0) error stop 'mixture response'
        call water%state_ph_checked(p+dp, h, a)
        call water%state_ph_checked(p-dp, h, b)
        if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'mixture pressure perturbation'
        call check_slope((a%properties%T-b%properties%T)/(2d0*dp), jac%ph(JAC_T,1))
        call check_slope((a%properties%rho-b%properties%rho)/(2d0*dp), jac%ph(JAC_RHO,1))
        call water%state_ps_checked(p+dp, s, a)
        call water%state_ps_checked(p-dp, s, b)
        if (a%status /= WATER_OK .or. b%status /= WATER_OK) error stop 'mixture entropy perturbation'
        call check_slope((a%properties%rho-b%properties%rho)/(2d0*dp), jac%ps(JAC_RHO,1))
    end subroutine verify_mixture

    subroutine check_slope(reference, actual)
        implicit none
        real(real64), intent(in) :: reference, actual
        if (abs(actual-reference) > 2d-14 + 2d-5*abs(reference)) then
            print *, reference, actual
            error stop 'mixture derivative'
        end if
    end subroutine check_slope
end program solver_jacobian
