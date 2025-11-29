program main
    use, intrinsic :: iso_fortran_env
    use :: module_iapws95
    implicit none

    type(type_iapws95) :: iapws95_model
    type(type_iapws95_phi0_properties) :: phi0_props
    type(type_iapws95_phir_properties) :: phir_props

    real(real64), parameter :: T = 500.0d0 ! Temperature in K
    real(real64), parameter :: rho = 838.025d0 ! Density in kg/m^3

    real(real64) :: tau, delta
    ! Calculate reduced properties
    tau = iapws95_model%T_c / T
    delta = rho / iapws95_model%rho_c

    ! Calculate ideal Helmholtz properties
    phi0_props = iapws95_model%calc_phi0_iapws95(tau, delta)
    ! Calculate residual Helmholtz properties
    phir_props = iapws95_model%calc_phir_iapws95(tau, delta)

    ! Output results
    print '(a10, es20.8)', "phi0: ", phi0_props%phi0
    print '(a10, es20.8)', "phi0_d: ", phi0_props%phi0_d
    print '(a10, es20.8)', "phi0_dd: ", phi0_props%phi0_dd
    print '(a10, es20.8)', "phi0_t: ", phi0_props%phi0_t
    print '(a10, es20.8)', "phi0_tt: ", phi0_props%phi0_tt
    print '(a10, es20.8)', "phi0_dt: ", phi0_props%phi0_dt
    print *, ""
    print '(a10, es20.8)', "phir: ", phir_props%phir
    print '(a10, es20.8)', "phir_d: ", phir_props%phir_d
    print '(a10, es20.8)', "phir_dd: ", phir_props%phir_dd
    print '(a10, es20.8)', "phir_t: ", phir_props%phir_t
    print '(a10, es20.8)', "phir_tt: ", phir_props%phir_tt
    print '(a10, es20.8)', "phir_dt", phir_props%phir_dt

end program main
