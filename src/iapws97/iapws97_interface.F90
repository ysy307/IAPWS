module module_iapws97
    use, intrinsic :: iso_fortran_env
    use :: module_iapws
    use :: utils_kahan, only:kahan_add
    implicit none
    private

    public :: type_iapws97

    type :: type_iapws97_auxiliary
        real(real64) :: T_star
        real(real64) :: p_star
        logical :: is_initialized = .false.
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_auxiliary
        procedure, pass(self), public :: calc_p_boundary => calc_p_boundary_iapws97_region23
        procedure, pass(self), public :: calc_T_boundary => calc_T_boundary_iapws97_region23
    end type type_iapws97_auxiliary

    interface
        module pure elemental subroutine initialize_type_iapws97_auxiliary(self)
            implicit none
            class(type_iapws97_auxiliary), intent(inout) :: self
        end subroutine initialize_type_iapws97_auxiliary

        module pure elemental function calc_p_boundary_iapws97_region23(self, temperature) result(pressure)
            implicit none
            class(type_iapws97_auxiliary), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64) :: pressure

        end function calc_p_boundary_iapws97_region23

        module pure elemental function calc_T_boundary_iapws97_region23(self, pressure) result(temperature)
            implicit none
            class(type_iapws97_auxiliary), intent(in) :: self
            real(real64), intent(in) :: pressure
            real(real64) :: temperature

        end function calc_T_boundary_iapws97_region23
    end interface

    type, extends(abst_iapws_gibbs) :: type_iapws97_region1
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_region1
        procedure, pass(self), public :: calc_gamma => calc_gamma_iapws97_region1
    end type type_iapws97_region1

    interface
        !> Initialize the IAPWS-97 Region 1 object.
        module pure elemental subroutine initialize_type_iapws97_region1(self)
            implicit none
            !> Initialize the IAPWS-97 Region 1 object.
            class(type_iapws97_region1), intent(inout) :: self
        end subroutine initialize_type_iapws97_region1

        !> Calculate the dimensionless Gibbs free energy \(\gamma\) for Region 1.
        module pure elemental subroutine calc_gamma_iapws97_region1(self, T_in, P_in, coef)
            implicit none
            class(type_iapws97_region1), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, p [Pa].
            real(real64), intent(in) :: P_in
            !> IAPWS Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        end subroutine calc_gamma_iapws97_region1
    end interface

    type, extends(abst_iapws_gibbs) :: type_iapws97_region2
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_region2
        procedure, pass(self), public :: calc_gamma => calc_gamma_iapws97_region2
    end type type_iapws97_region2

    interface
        !> Initialize the IAPWS-97 Region 2 object.
        module pure elemental subroutine initialize_type_iapws97_region2(self)
            implicit none
            !> IAPWS-97 Region 2 object.
            class(type_iapws97_region2), intent(inout) :: self
        end subroutine initialize_type_iapws97_region2

        !> Calculate the Gibbs free energy \(\gamma\) for Region 2.
        module pure elemental subroutine calc_gamma_iapws97_region2(self, T_in, P_in, coef)
            implicit none
            !> IAPWS-97 Region 2 object.
            class(type_iapws97_region2), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, p [Pa].
            real(real64), intent(in) :: P_in
            !> IAPWS Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        end subroutine calc_gamma_iapws97_region2
    end interface

    type, extends(abst_iapws_helmholtz) :: type_iapws97_region3
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_region3
        procedure, pass(self), public :: calc_phi => calc_phi_iapws97_region3
    end type type_iapws97_region3

    interface
        module pure elemental subroutine initialize_type_iapws97_region3(self)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws97_region3), intent(inout) :: self
        end subroutine initialize_type_iapws97_region3

        module pure elemental subroutine calc_phi_iapws97_region3(self, tau, delta, property)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws97_region3), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 helmholtz properties
            type(type_iapws_helmholtz_property), intent(inout) :: property
        end subroutine calc_phi_iapws97_region3
    end interface

    type :: type_iapws97_region4
        real(real64) :: T_star
        real(real64) :: p_star
        logical :: is_initialized = .false.
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_region4
        procedure, pass(self), public :: calc_psat => calc_psat_iapws97_region4
        procedure, pass(self), public :: calc_tsat => calc_tsat_iapws97_region4
    end type type_iapws97_region4

    interface
        module pure elemental subroutine initialize_type_iapws97_region4(self)
            implicit none
            class(type_iapws97_region4), intent(inout) :: self

        end subroutine initialize_type_iapws97_region4

        module pure elemental function calc_psat_iapws97_region4(self, temperature) result(P_sat)
            implicit none
            class(type_iapws97_region4), intent(in) :: self
            real(real64), intent(in) :: temperature
            real(real64) :: P_sat

        end function calc_psat_iapws97_region4

        module pure elemental function calc_tsat_iapws97_region4(self, pressure) result(T_sat)
            implicit none
            class(type_iapws97_region4), intent(in) :: self
            real(real64), intent(in) :: pressure
            real(real64) :: T_sat

        end function calc_tsat_iapws97_region4
    end interface

    type, extends(abst_iapws_gibbs) :: type_iapws97_region5
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97_region5
        procedure, pass(self), public :: calc_gamma => calc_gamma_iapws97_region5
    end type type_iapws97_region5

    interface
        !> Initialize the IAPWS-97 Region 5 object.
        module pure elemental subroutine initialize_type_iapws97_region5(self)
            implicit none
            !> IAPWS-97 Region 5 object.
            class(type_iapws97_region5), intent(inout) :: self
        end subroutine initialize_type_iapws97_region5

        !> Calculate the dimensionless Gibbs free energy \(\gamma\) for Region 5.
        module pure elemental subroutine calc_gamma_iapws97_region5(self, T_in, P_in, coef)
            implicit none
            !> IAPWS-97 Region 5 object.
            class(type_iapws97_region5), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, p [Pa].
            real(real64), intent(in) :: P_in
            !> IAPWS Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        end subroutine calc_gamma_iapws97_region5
    end interface

    type :: type_iapws97
        private
        type(type_iapws97_auxiliary) :: auxiliary
        type(type_iapws97_region1) :: region1
        type(type_iapws97_region2) :: region2
        type(type_iapws97_region3) :: region3
        type(type_iapws97_region4) :: region4
        type(type_iapws97_region5) :: region5
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws97
        procedure, pass(self), public :: get_region => get_region_iapws97
        procedure, pass(self), public :: calc_properties => calc_properties_iapws97
        procedure, pass(self), public :: calc_nu => calc_nu_iapws97
        procedure, pass(self), public :: calc_rho => calc_rho_iapws97
        procedure, pass(self), public :: calc_drho_dT => calc_drho_dT_iapws97
        procedure, pass(self), public :: calc_drho_dp => calc_drho_dp_iapws97
        procedure, pass(self), public :: calc_u => calc_u_iapws97
        procedure, pass(self), public :: calc_h => calc_h_iapws97
        procedure, pass(self), public :: calc_s => calc_s_iapws97
        procedure, pass(self), public :: calc_cp => calc_cp_iapws97
        procedure, pass(self), public :: calc_cv => calc_cv_iapws97
        procedure, pass(self), public :: calc_w => calc_w_iapws97
        procedure, pass(self), public :: calc_latent_heat => calc_latent_heat_iapws97
    end type type_iapws97

    interface
        module pure elemental subroutine initialize_type_iapws97(self)
            implicit none
            class(type_iapws97), intent(inout) :: self
        end subroutine initialize_type_iapws97

        module pure elemental function get_region_iapws97(self, T_in, p_in) result(region_id)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            integer(int32) :: region_id

        end function get_region_iapws97

        module pure elemental subroutine calc_properties_iapws97(self, T_in, p_in, property)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            type(type_iapws_property), intent(inout) :: property

        end subroutine calc_properties_iapws97

        module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: nu

        end subroutine calc_nu_iapws97

        module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: rho

        end subroutine calc_rho_iapws97

        module pure elemental subroutine calc_drho_dT_iapws97(self, T_in, p_in, drho_dT)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: drho_dT

        end subroutine calc_drho_dT_iapws97

        module pure elemental subroutine calc_drho_dp_iapws97(self, T_in, p_in, drho_dp)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: drho_dp

        end subroutine calc_drho_dp_iapws97

        module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: u

        end subroutine calc_u_iapws97

        module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: h

        end subroutine calc_h_iapws97

        module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: s

        end subroutine calc_s_iapws97

        module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: cp

        end subroutine calc_cp_iapws97

        module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: cv

        end subroutine calc_cv_iapws97

        module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: w

        end subroutine calc_w_iapws97

        !> Calculate Latent Heat of Vaporization [J/kg]
        !> User must provide EITHER T_in [K] OR p_in [Pa].
        module pure elemental subroutine calc_latent_heat_iapws97(self, latent_heat, T_in, p_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(inout) :: latent_heat
            real(real64), intent(in), optional :: T_in
            real(real64), intent(in), optional :: p_in
        end subroutine calc_latent_heat_iapws97
    end interface

end module module_iapws97
