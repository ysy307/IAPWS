module module_iapws95
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws, only:abst_iapws_helmholtz, type_iapws_phi_property
    use :: iapws95_constants
    implicit none
    private

    public :: type_iapws95
    ! public :: type_iapws95_phi0_properties
    ! public :: type_iapws95_phir_properties

    type, extends(abst_iapws_helmholtz) :: type_iapws95
        ! real(real64) :: T_c = critical_temperature
        ! real(real64) :: rho_c = critical_density
    contains
        procedure, pass(self) :: initialize => initialize_type_iapws95
        procedure, pass(self) :: calc_phi => calc_phi_iapws95
        procedure, pass(self) :: calc_phi0 => calc_phi0_iapws95
        procedure, pass(self) :: calc_phir => calc_phir_iapws95
    end type type_iapws95

    interface
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
            type(type_iapws_phi_property), intent(inout) :: property
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
            type(type_iapws_phi_property), intent(inout) :: property

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
            type(type_iapws_phi_property), intent(inout) :: property

        end subroutine calc_phir_iapws95
    end interface

end module module_iapws95
