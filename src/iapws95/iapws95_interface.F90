module module_iapws95
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use :: iapws95_constants
    implicit none

    public :: type_iapws95
    public :: type_iapws95_phi0_properties
    public :: type_iapws95_phir_properties

    type :: type_iapws95
        real(real64) :: T_c = critical_temperature
        real(real64) :: rho_c = critical_density
    contains
        procedure, nopass, public :: calc_phi0_iapws95
        procedure, nopass, public :: calc_phir_iapws95
    end type type_iapws95

    !> Data structure to hold IAPWS-95 ideal properties
    type :: type_iapws95_phi0_properties
        !> Ideal adimensional helmholtz energy [-]
        real(real64) :: phi0
        !> Ideal adimensional helmholtz energy derivative with respect to delta,∂phi0/∂δ|τ [-]
        real(real64) :: phi0_d
        !> Ideal adimensional helmholtz energy derivative with respect to tau,∂phi0/∂τ|δ [-]
        real(real64) :: phi0_t
        !> Ideal adimensional helmholtz energy second derivative with respect to delta,∂²phi0/∂δ²|τ [-]
        real(real64) :: phi0_dd
        !> Ideal adimensional helmholtz energy second derivative with respect to tau,∂²phi0/∂τ²|δ [-]
        real(real64) :: phi0_tt
        !> Ideal adimensional helmholtz energy mixed second derivative,∂²phi0/∂δ∂τ [-]
        real(real64) :: phi0_dt
    end type type_iapws95_phi0_properties

    !> Data structure to hold IAPWS-95 ideal properties
    type :: type_iapws95_phir_properties
        !> Ideal adimensional helmholtz energy [-]
        real(real64) :: phir
        !> Ideal adimensional helmholtz energy derivative with respect to delta,∂phir/∂δ|τ [-]
        real(real64) :: phir_d
        !> Ideal adimensional helmholtz energy derivative with respect to tau,∂phir/∂τ|δ [-]
        real(real64) :: phir_t
        !> Ideal adimensional helmholtz energy second derivative with respect to delta,∂²phir/∂δ²|τ [-]
        real(real64) :: phir_dd
        !> Ideal adimensional helmholtz energy second derivative with respect to tau,∂²phir/∂τ²|δ [-]
        real(real64) :: phir_tt
        !> Ideal adimensional helmholtz energy mixed second derivative,∂²phir/∂δ∂τ [-]
        real(real64) :: phir_dt
    end type type_iapws95_phir_properties

    interface
        module pure elemental function calc_phi0_iapws95(tau, delta) result(property)
            implicit none
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws95_phi0_properties) :: property

        end function calc_phi0_iapws95

        module pure elemental function calc_phir_iapws95(tau, delta) result(property)
            implicit none
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws95_phir_properties) :: property

        end function calc_phir_iapws95
    end interface

contains
    function example_function(x) result(y)
        real(real64), intent(in) :: x
        real(real64) :: y

        ! Example implementation: y = x^2
        y = x * x
    end function example_function

end module module_iapws95
