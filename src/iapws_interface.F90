module iapws
    use, intrinsic :: iso_fortran_env
    implicit none

    type, abstract :: abst_iapws_helmholtz
        !> Critical Temperature, \( T_c \) [K]
        real(real64) :: T_c
        !> Critical Density, \( \rho_c \) [kg/m^3]
        real(real64) :: rho_c
        !> Specific Gas Constant, \( R \) [J/(kg K)]
        real(real64) :: R
        !> Initialization flag
        logical :: is_initialized = .false.
    contains
        procedure(abst_iapws_helmholtz_initialize), pass(self), public, deferred :: initialize !&
        procedure(abst_calc_phi_iapws),             pass(self), public, deferred :: calc_phi !&
        procedure, pass(self), public :: calc_properties => calc_properties_helmholtz
        procedure, pass(self), public :: calc_p => calc_p_helmholtz
        procedure, pass(self), public :: calc_u => calc_u_helmholtz
        procedure, pass(self), public :: calc_s => calc_s_helmholtz
        procedure, pass(self), public :: calc_h => calc_h_helmholtz
        procedure, pass(self), public :: calc_cp => calc_cp_helmholtz
        procedure, pass(self), public :: calc_cv => calc_cv_helmholtz
        procedure, pass(self), public :: calc_w => calc_w_helmholtz
    end type abst_iapws_helmholtz

    type, abstract :: abst_iapws_gibbs_model
        integer(int32) :: T_c
        integer(int32) :: rho_c
    contains
        ! Gibbs related procedures can be added here in the future
    end type abst_iapws_gibbs_model

    !> Data structure to hold IAPWS properties
    type :: type_iapws_phi_property
        !> Total adimensional helmholtz energy [-]
        real(real64) :: phi
        !> Derivative of total adimensional helmholtz energy with respect to delta,∂phi/∂δ|τ [-]
        real(real64) :: phi_d
        !> Derivative of total adimensional helmholtz energy with respect to tau,∂phi/∂τ|δ [-]
        real(real64) :: phi_t
        !> Second derivative of total adimensional helmholtz energy with respect to delta,∂²phi/∂δ²|τ [-]
        real(real64) :: phi_dd
        !> Second derivative of total adimensional helmholtz energy with respect to tau,∂²phi/∂τ²|δ [-]
        real(real64) :: phi_tt
        !> Mixed second derivative of total adimensional helmholtz energy,∂²phi/∂δ∂τ [-]
        real(real64) :: phi_dt
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
    end type type_iapws_phi_property

    !> IAPWS properties data structure
    type :: type_iapws_property
        !> Region ID
        integer(int32) :: region_id = -1
        !> Specific Volume, \( \nu \) [m^3/kg]
        real(real64) :: nu = 0.0d0
        !> Density, \( \rho \) [kg/m^3]
        real(real64) :: rho = 0.0d0
        !> Specific Internal Energy, \( u \) [J/kg]
        real(real64) :: u = 0.0d0
        !> Specific Entropy, \( s \) [J/(kg K)]
        real(real64) :: s = 0.0d0
        !> Specific Enthalpy, \( h \) [J/kg]
        real(real64) :: h = 0.0d0
        !> Specific Heat Capacity at constant pressure, \( c_p \) [J/(kg K)]
        real(real64) :: cp = 0.0d0
        !> Specific Heat Capacity at constant volume, \( c_v \) [J/(kg K)]
        real(real64) :: cv = 0.0d0
        !> Speed of Sound, \( w \) [m/s]
        real(real64) :: w = 0.0d0
        !> Pressure, \( p \) [Pa]
        real(real64) :: p = 0.0d0
        !> Temperature, \( T \) [K]
        real(real64) :: T = 0.0d0
        !> Cubic expansion coefficient [1/K]
        real(real64) :: alpha = 0.0d0
        !> Pressure coefficient [Pa/K]
        real(real64) :: beta = 0.0d0
        !> Isothermal compressibility [1/Pa]
        real(real64) :: kappa_s = 0.0d0
        !> Isentropic compressibility [1/Pa]
        real(real64) :: kappa_T = 0.0d0
    end type type_iapws_property

    abstract interface
        !> Initialize IAPWS Helmholtz model
        pure elemental subroutine abst_iapws_helmholtz_initialize(self)
            import :: abst_iapws_helmholtz
            implicit none
            !> IAPWS Helmholtz model instance
            class(abst_iapws_helmholtz), intent(inout) :: self
        end subroutine abst_iapws_helmholtz_initialize

        !> Calculate IAPWS Helmholtz properties
        pure elemental subroutine abst_calc_phi_iapws(self, tau, delta, property)
            import :: abst_iapws_helmholtz, type_iapws_phi_property, real64
            implicit none
            !> IAPWS Helmholtz model instance
            class(abst_iapws_helmholtz), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws_phi_property), intent(inout) :: property
        end subroutine abst_calc_phi_iapws
    end interface

    interface
        module pure elemental subroutine calc_properties_helmholtz(self, T_in, rho_in, property)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            type(type_iapws_property), intent(inout) :: property

        end subroutine calc_properties_helmholtz

        module pure elemental subroutine calc_p_helmholtz(self, T_in, rho_in, p, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: p
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_p_helmholtz

        module pure elemental subroutine calc_u_helmholtz(self, T_in, rho_in, u, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: u
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_u_helmholtz

        module pure elemental subroutine calc_s_helmholtz(self, T_in, rho_in, s, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: s
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_s_helmholtz

        module pure elemental subroutine calc_h_helmholtz(self, T_in, rho_in, h, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: h
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_h_helmholtz

        module pure elemental subroutine calc_cp_helmholtz(self, T_in, rho_in, cp, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: cp
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_cp_helmholtz

        module pure elemental subroutine calc_cv_helmholtz(self, T_in, rho_in, cv, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: cv
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_cv_helmholtz

        module pure elemental subroutine calc_w_helmholtz(self, T_in, rho_in, w, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: w
            type(type_iapws_phi_property), intent(inout), optional :: prop_in

        end subroutine calc_w_helmholtz
    end interface

end module iapws
