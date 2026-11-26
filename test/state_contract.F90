program state_contract
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iapws
    implicit none
    type(type_iapws95) :: scientific
    type(type_iapws97) :: industrial
    type(type_iapws_property) :: state
    real(real64) :: rho, p, batch_rho(1)
    integer :: status(1)

    call scientific%initialize()
    call scientific%calc_rho(300.0_real64, 1.0e5_real64, rho)
    call scientific%calc_p(300.0_real64, rho, p)
    if (abs(rho - 996.5563403889_real64) > 1.0e-8_real64) error stop 'IAPWS95 liquid density inverse'
    if (abs(p - 1.0e5_real64) > 1.0e-5_real64) error stop 'IAPWS95 inverse pressure residual'

    call industrial%initialize()
    call industrial%calc_properties(300.0_real64, 1.0e5_real64, state)
    if (state%T /= 300.0_real64 .or. state%p /= 1.0e5_real64) error stop 'Gibbs state coordinates'
    call industrial%calc_properties(650.0_real64, 25.0e6_real64, state)
    if (abs(state%nu * state%rho - 1.0_real64) > epsilon(rho)) error stop 'Helmholtz specific volume'
    if (abs(state%alpha + state%drho_dT / state%rho) > epsilon(rho)) error stop 'Helmholtz expansivity'

    call industrial%calc_properties_batch([647.095990_real64], [22063997.3189_real64], &
                                          rho=batch_rho, status=status)
    if (ieee_is_finite(batch_rho(1)) .neqv. (status(1) == IAPWS_BATCH_OK)) error stop 'IF97 batch solve status'
end program state_contract
