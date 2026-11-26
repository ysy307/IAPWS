program water_state
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use iapws
    implicit none
    type(type_water) :: water
    type(type_iapws97) :: industrial
    type(type_iapws95) :: scientific
    type(type_iapws_property) :: exact_state
    type(type_water_state) :: state, state2
    type(type_water_saturation) :: sat, sat_p
    real(real64), parameter :: t(5) = [260d0, 254d0, 265d0, 320d0, 550d0]
    ! iapws 1.5.5 scipy/numpy reference evaluation, in Pa.
    real(real64), parameter :: p(5) = [138268113.00221789d0, 268684646.63361079d0, &
        479640244.37879908d0, 1356756517.8693883d0, 6308714243.5430183d0]
    ! R14-08(2011) table 3, in Pa.
    real(real64), parameter :: rounded(5) = [138.268d6, 268.685d6, 479.640d6, 1356.76d6, 6308.71d6]
    integer, parameter :: solids(5) = [ICE_IH, ICE_III, ICE_V, ICE_VI, ICE_VII]
    integer :: i
    real(real64) :: pm, mu, conductivity, slope

    do i = 1, 5
        pm = melting_pressure(t(i), solids(i))
        if (abs(pm / p(i) - 1d0) > 16d0 * epsilon(pm)) error stop 'Python melting pressure'
        if (abs(pm - rounded(i)) > merge(500d0, 5000d0, i <= 3)) error stop 'R14 table 3 melt'
    end do
    pm = sublimation_pressure(230d0)
    if (abs(pm / 8.9473527401891513d0 - 1d0) > 16d0 * epsilon(pm)) error stop 'Python sublimation pressure'
    if (abs(pm - 8.94735d0) > 5d-6) error stop 'R14 table 3 sublimation'

    water = new_water(fluid=IF97)
    call water%state_tp(300d0, 1d5, state)
    if (state%status /= WATER_OK .or. state%phase /= LIQUID .or. state%formulation /= IF97) error stop 'IF97 liquid'
    call water%state_tp(260d0, 1d5, state)
    if (state%status /= WATER_OK .or. state%phase /= ICE_IH .or. state%formulation /= IAPWS06) error stop 'IF97 ice'
    call water%state_tp(260d0, p(1), state)
    if (state%status /= WATER_COEXISTENCE .or. .not. ieee_is_nan(state%properties%rho)) error stop 'melting coexistence'
    call water%state_tp(260d0, 150d6, state)
    if (state%status /= WATER_OUT_OF_RANGE .or. state%phase /= LIQUID) error stop 'cold stable liquid outside IF97'
    water = new_water(fluid=IF97, extend_cold_liquid=.true.)
    call water%state_tp(260d0, 150d6, state)
    if (state%status /= WATER_OK .or. state%formulation /= IAPWS95) error stop 'explicit cold-liquid extension'
    if (abs(state%properties%rho - 1067.9400725589032d0) > 2d-9) error stop 'cold density Python reference'
    call water%state_tp(254d0, 300d6, state)
    if (state%status /= WATER_UNSUPPORTED .or. state%phase /= HIGH_PRESSURE_ICE) error stop 'high pressure ice'
    water = new_water(fluid=IAPWS95)
    call water%saturation_t(450d0, sat)
    if (sat%status /= WATER_OK) error stop 'IAPWS95 saturation'
    if (abs(sat%p - 932203.5636281562d0) > 2d-7) error stop 'IAPWS95 saturation pressure'
    if (abs(sat%liquid%rho - 890.3412497616729d0) > 2d-10) error stop 'IAPWS95 saturation liquid density'
    call water%state_tp(235.15d0, 101325d0, state, phase=LIQUID, allow_metastable=.true.)
    if (state%status /= WATER_METASTABLE .or. state%formulation /= G12_15) error stop 'G12 metastable selection'
    if (abs(state%properties%rho - 968.0999868967568d0) > 1d-9) error stop 'G12 Python density'
    if (abs(state%properties%cp - 5997.563d0) > 5d-4) error stop 'G12 table 5 heat capacity'
    call water%state_tp(230d0, 101325d0, state, phase=LIQUID, allow_metastable=.true.)
    if (state%status /= WATER_OUT_OF_RANGE) error stop 'G12 nucleation limit'
    call water%state_tp(100d0, 1d-16, state)
    if (state%status /= WATER_OK .or. state%formulation /= G9_12) error stop 'G9 low temperature vapor'
    if (abs(state%properties%cp / 461.51805d0 - 4.00536708d0) > 5d-9) error stop 'G9 table 3 ideal heat capacity'
    water = new_water(fluid=IF97)
    call water%state_tp(700d0, 3d6, state)
    call water%state_ph(3d6, state%properties%h, state2)
    if (state2%status /= WATER_OK .or. abs(state2%properties%T - 700d0) > 4d-12) error stop 'IF97 p,h round trip'
    call water%state_ps(3d6, state%properties%s, state2)
    if (state2%status /= WATER_OK .or. abs(state2%properties%T - 700d0) > 4d-12) error stop 'IF97 p,s round trip'
    call water%saturation_t(500d0, sat)
    call water%state_ph(sat%p, (sat%liquid%h + sat%vapor%h) * 0.5d0, state)
    if (state%status /= WATER_OK .or. abs(state%quality - 0.5d0) > 2d-14) error stop 'IF97 wet steam p,h'
    mu = viscosity(298.15d0, 998d0)
    conductivity = thermal_conductivity(298.15d0, 998d0)
    if (abs(mu / 0.0008897351001498108d0 - 1d0) > 128d0 * epsilon(mu)) error stop 'R12 Python viscosity'
    if (abs(mu * 1d6 - 889.735100d0) > 5d-7) error stop 'R12 table 4 viscosity'
    if (abs(conductivity / 0.6077128675880629d0 - 1d0) > 128d0 * epsilon(mu)) error stop 'R15 Python conductivity'
    if (abs(conductivity * 1d3 - 607.712868d0) > 5d-7) error stop 'R15 table 4 conductivity'
    call scientific%initialize()
    call scientific%calc_properties(647.35d0, 322d0, exact_state)
    call scientific%calc_drho_dp(970.644d0, 322d0, slope)
    mu = viscosity(647.35d0, 322d0, exact_state%drho_dP, slope)
    conductivity = thermal_conductivity(647.35d0, 322d0, exact_state%cp, exact_state%cv, mu, exact_state%drho_dP, slope)
    if (abs(mu * 1d6 - 42.961579d0) > 5d-7) error stop 'R12 table 5 critical viscosity'
    if (abs(conductivity * 1d3 - 1443.75556d0) > 5d-6) error stop 'R15 table 5 critical conductivity'
    if (abs(surface_tension(300d0) / 0.07168596252716256d0 - 1d0) > 16d0 * epsilon(mu)) &
        error stop 'R1 Python surface tension'
    water = new_water(fluid=IAPWS95)
    call water%saturation_t(450d0, sat)
    call water%state_tp(450d0, sat%p, state)
    if (state%status /= WATER_COEXISTENCE) error stop '95 saturation coexistence'
    call water%state_tp(450d0, sat%p, state, phase=VAPOR)
    if (state%status /= WATER_OK .or. state%properties%rho /= sat%vapor%rho) error stop '95 saturation phase'
    call water%state_tp(260d0, 1d5, state)
    call water%state_ph(1d5, state%properties%h, state2)
    if (state2%status /= WATER_OK .or. state2%phase /= ICE_IH) error stop '95 ice p,h'
    if (abs(state2%properties%T - 260d0) > 4d-12) error stop '95 ice temperature'
    call water%state_ph(1d5, -1.5d5, state2)
    if (state2%status /= WATER_OK .or. state2%phase /= ICE_LIQUID) error stop '95 melting p,h'
    if (abs(state2%properties%T - 273.15261752072672d0) > 5d-7) error stop '95 melting temperature'
    if (state2%solid_fraction <= 0d0 .or. state2%solid_fraction >= 1d0) error stop '95 melting fraction'
    call industrial%initialize()
    if (abs(industrial%backward_temperature(3d6, 500d3, 1) - 391.79850876242557d0) > 2d-12) &
        error stop 'IF97 backward T(p,h) region 1'
    if (abs(industrial%backward_temperature(80d6, 3d3, 1, entropy=.true.) - 565.8999087585007d0) > 2d-12) &
        error stop 'IF97 backward T(p,s) region 1'
    if (abs(industrial%backward_temperature(5d6, 3500d3, 2) - 801.299101856204d0) > 3d-12) &
        error stop 'IF97 backward T(p,h) region 2'
    if (abs(industrial%backward_temperature(0.1d6, 7.5d3, 2, entropy=.true.) - 399.51709665541114d0) > 3d-12) &
        error stop 'IF97 backward T(p,s) region 2'
    if (abs(industrial%backward_temperature(25d6, 2000d3, 3) - 655.3573690102869d0) > 3d-12) &
        error stop 'IF97 backward T(p,h) region 3'
    if (abs(industrial%backward_temperature(25d6, 4d3, 3, entropy=.true.) - 646.4355837040306d0) > 3d-12) &
        error stop 'IF97 backward T(p,s) region 3'
    call water%state_tp(300d0, 1d5, state)
    call water%transport(state, mu=mu, conductivity=conductivity)
    if (abs(mu / 0.0008537426062537289d0 - 1d0) > 1d-11) error stop '95 state viscosity'
    if (abs(conductivity / 0.6094991151038355d0 - 1d0) > 1d-11) error stop '95 state conductivity'
    water = new_water(fluid=IF97)
    call water%state_tp(450d0, 1d6, state, phase=VAPOR, allow_metastable=.true.)
    if (state%status /= WATER_METASTABLE .or. state%formulation /= IF97_METASTABLE_VAPOR) error stop '97 metastable vapor'
    if (abs(state%properties%nu / 0.1925165400597291d0 - 1d0) > 16d0 * epsilon(mu)) error stop '97 meta Python'
    if (abs(state%properties%nu - 0.192516540d0) > 5d-10) error stop '97 meta Table 18'
    call water%state_tp(350d0, 1d6, state, phase=VAPOR, allow_metastable=.true.)
    if (state%status /= WATER_OUT_OF_RANGE) error stop '97 metastable moisture limit'
    call water%state_tp(273.16d0, 611.657d0, state, phase=LIQUID)
    if (state%status /= WATER_OK .or. state%properties%rho < 900d0) error stop 'IF97 triple liquid'
    water = new_water(fluid=IAPWS95)
    call water%state_tp(273.16d0, 611.657d0, state, phase=VAPOR)
    if (state%status /= WATER_OK .or. state%properties%rho > 0.01d0) error stop '95 triple vapor'
    call water%saturation_t(273.16d0, sat)
    call water%state_tp(273.16d0, sat%p, state)
    if (state%status /= WATER_COEXISTENCE) error stop '95 numerical triple saturation'
    do i = 1, 2
        water = new_water(fluid=merge(IF97, IAPWS95, i == 1))
        call water%state_tp(300d0, 22063999d0, state)
        call water%state_ph(22063999d0, state%properties%h, state2)
        if (state2%status /= WATER_OK) error stop 'near-critical isobar away from saturation'
        if (abs(state2%properties%T - 300d0) > 2d-9) error stop 'near-critical isobar inversion'
    end do
    water = new_water(fluid=IF97)
    call water%state_tp(1073.15d0, 25d6, state)
    call water%state_ps(25d6, state%properties%s, state2)
    if (state2%status /= WATER_OK .or. state2%properties%T /= 1073.15d0) error stop 'IF97 Region 2/5 seam'
    call water%state_tp(623.15d0, 25d6, state)
    call water%state_ph(25d6, state%properties%h, state2)
    if (state2%status /= WATER_OK .or. state2%properties%T /= 623.15d0) error stop 'IF97 Region 1/3 seam'
    call water%state_tp(647.096d0, 22.064d6, state)
    if (state%status /= WATER_SINGULAR .or. state%properties%rho /= 322d0) error stop 'critical TP'
    if (.not. ieee_is_nan(state%properties%cp)) error stop 'critical cp'
    water = new_water(fluid=IAPWS95)
    call water%state_tp(254d0, p(2), state)
    if (state%status /= WATER_COEXISTENCE) error stop 'high pressure melting coexistence'
    call water%state_tp(254d0, p(2), state, phase=LIQUID)
    if (state%status /= WATER_OK .or. state%formulation /= IAPWS95) error stop 'high pressure melting liquid'
    call water%state_tp(300d0, 611.757d0, state)
    call water%state_ps(611.757d0, state%properties%s, state2)
    if (state2%status /= WATER_OK) error stop 'saturation inverse at low pressure'
    call water%saturation_t(500d0, sat)
    call water%saturation_p(sat%p, sat_p)
    if (sat_p%status /= WATER_OK .or. abs(sat_p%T - 500d0) > 2d-10) error stop '95 saturation p-to-T'
    water = new_water(fluid=IF97)
    call water%saturation_t(500d0, sat)
    call water%saturation_p(sat%p, sat_p)
    if (sat_p%status /= WATER_OK .or. abs(sat_p%T - 500d0) > 2d-10) error stop 'IF97 saturation p-to-T'
    call water%state_tp(1073.17d0, 25d6, state)
    call water%state_ps(25d6, state%properties%s, state2)
    if (state2%status /= WATER_AMBIGUOUS) error stop 'IF97 seam overlap'
    call water%state_ps(25d6, state%properties%s, state2, region=5)
    if (state2%status /= WATER_OK .or. abs(state2%properties%T - 1073.17d0) > 5d-10) error stop 'IF97 region hint'
    water = new_water(fluid=IAPWS95)
    call water%saturation_t(390.98417495761834d0, sat)
    if (sat%status /= WATER_OK) error stop '95 saturation inside ordinary range'
    ! iapws 1.5.5 coexistence at the temperature that previously failed.
    if (abs(sat%p / 185422.64322053714d0 - 1d0) > 8d-12) error stop '95 ordinary saturation pressure'
    if (abs(sat%liquid%rho / 944.8417464982938d0 - 1d0) > 8d-12) error stop '95 ordinary saturation liquid'
    if (abs(sat%vapor%rho / 1.0518618480753084d0 - 1d0) > 8d-12) error stop '95 ordinary saturation vapor'
end program water_state
