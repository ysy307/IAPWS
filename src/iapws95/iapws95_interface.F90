module module_iapws95
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: module_iapws, only:abst_iapws_helmholtz, type_iapws_helmholtz_property, type_iapws_property
    use :: iapws95_constants
    implicit none
    private

    public :: type_iapws95

    type, extends(abst_iapws_helmholtz) :: type_iapws95
    contains
        procedure, pass(self) :: initialize => initialize_type_iapws95
        procedure, pass(self) :: saturation_estimate => saturation_estimate_iapws95
        procedure, pass(self) :: saturation_densities => saturation_densities_iapws95
        procedure, pass(self) :: solve_density => solve_density_iapws95
        procedure, pass(self) :: calc_rho => calc_rho_iapws95
        procedure, pass(self) :: calc_phi => calc_phi_iapws95
        procedure, pass(self) :: calc_phi0 => calc_phi0_iapws95
        procedure, pass(self) :: calc_phir => calc_phir_iapws95
    end type type_iapws95

    interface
        module pure elemental function saturation_estimate_iapws95(self, T) result(p)
            implicit none
            real(real64) :: p
            class(type_iapws95), intent(in) :: self
            real(real64), intent(in) :: T
        end function saturation_estimate_iapws95

        module pure elemental subroutine saturation_densities_iapws95(self, T, p, liquid, vapor, liquid_state, vapor_state)
            implicit none
            class(type_iapws95), intent(in) :: self
            real(real64), intent(in) :: T
            real(real64), intent(out) :: p, liquid, vapor
            type(type_iapws_property), intent(inout), optional :: liquid_state, vapor_state
        end subroutine saturation_densities_iapws95

        module pure elemental subroutine solve_density_iapws95(self, T_in, P_in, rho, guess_rho)
            implicit none
            class(type_iapws95), intent(in) :: self
            real(real64), intent(in) :: T_in, P_in
            real(real64), intent(out) :: rho
            real(real64), intent(in), optional :: guess_rho
        end subroutine solve_density_iapws95

        module pure elemental subroutine calc_rho_iapws95(self, T_in, P_in, rho)
            implicit none
            class(type_iapws95), intent(in) :: self
            real(real64), intent(in) :: T_in, P_in
            real(real64), intent(inout) :: rho
        end subroutine calc_rho_iapws95

        module pure elemental subroutine initialize_type_iapws95(self)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws95), intent(inout) :: self
        end subroutine initialize_type_iapws95

        module pure elemental subroutine calc_phi_iapws95(self, tau, delta, property)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws95), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 helmholtz properties
            type(type_iapws_helmholtz_property), intent(inout) :: property
        end subroutine calc_phi_iapws95

        module pure elemental subroutine calc_phi0_iapws95(self, tau, delta, property)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws95), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws_helmholtz_property), intent(inout) :: property

        end subroutine calc_phi0_iapws95

        module pure elemental subroutine calc_phir_iapws95(self, tau, delta, property)
            implicit none
            !> IAPWS-95 model instance
            class(type_iapws95), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws_helmholtz_property), intent(inout) :: property

        end subroutine calc_phir_iapws95
    end interface

end module module_iapws95
