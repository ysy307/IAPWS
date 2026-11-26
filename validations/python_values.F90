! Emit SI values for an independent Python calculation of the same states.
program python_values
    use, intrinsic :: iso_fortran_env, only: real64, iostat_end
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use iapws
    implicit none
    type(type_water) :: scientific, industrial
    type(type_water_state) :: state
    type(type_water_saturation) :: sat
    character(10) :: model
    character(25) :: mu_text, k_text
    real(real64) :: T, p, mu, k
    integer :: ios
    scientific = new_water(IAPWS95)
    industrial = new_water(IF97)
    do
        read (*,*,iostat=ios) model, T, p
        if (ios == iostat_end) exit
        if (ios /= 0) error stop 'invalid Python verification request'
        mu = 0d0
        k = 0d0
        select case (model)
        case ('95')
            call scientific%state_tp(T, p, state)
            call scientific%transport(state, mu, k)
        case ('97')
            call industrial%state_tp(T, p, state)
            call industrial%transport(state, mu, k)
        case ('ice')
            call industrial%state_tp(T, p, state)
        case ('meta')
            call industrial%state_tp(T, p, state, phase=VAPOR, allow_metastable=.true.)
        case ('g12')
            call scientific%state_tp(T, p, state, phase=LIQUID, allow_metastable=.true.)
        case ('g9')
            call scientific%state_tp(T, p, state)
        case ('sat')
            call scientific%saturation_t(T, sat)
            if (sat%status /= WATER_OK) error stop 'invalid verification saturation state'
            write (*,'(*(es25.16,1x))') sat%p, sat%liquid%rho, sat%vapor%rho, &
                                       sat%liquid%h, sat%vapor%h, sat%liquid%s, sat%vapor%s
            cycle
        case default
            error stop 'unknown Python verification model'
        end select
        if (state%status /= WATER_OK .and. state%status /= WATER_METASTABLE) then
            write (*,*) model, T, p, state%status, state%phase
            error stop 'invalid verification state'
        end if
        mu_text = 'nan'
        k_text = 'nan'
        if (.not. ieee_is_nan(mu)) write (mu_text,'(es25.16)') mu
        if (.not. ieee_is_nan(k)) write (k_text,'(es25.16)') k
        write (*,'(5(es25.16,1x),2(a,1x),2(es25.16,1x))') state%properties%rho, state%properties%h, state%properties%s, &
            state%properties%cp, state%properties%w, trim(mu_text), trim(k_text), state%properties%alpha, state%properties%kappa_T
    end do
end program python_values
