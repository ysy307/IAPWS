module module_iapws
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
        procedure, pass(self), public :: calc_rho => calc_rho_helmholtz
        procedure, pass(self), public :: calc_p => calc_p_helmholtz
        procedure, pass(self), public :: calc_p_rho => calc_p_rho_helmholtz
        procedure, pass(self), public :: calc_p_T => calc_p_T_helmholtz
        procedure, pass(self), public :: calc_u => calc_u_helmholtz
        procedure, pass(self), public :: calc_s => calc_s_helmholtz
        procedure, pass(self), public :: calc_h => calc_h_helmholtz
        procedure, pass(self), public :: calc_cp => calc_cp_helmholtz
        procedure, pass(self), public :: calc_cv => calc_cv_helmholtz
        procedure, pass(self), public :: calc_w => calc_w_helmholtz
    end type abst_iapws_helmholtz

    type, abstract :: abst_iapws_gibbs
        !> Reference Temperature, \( T_{star} \) [K]
        real(real64) :: T_star
        !> Reference Pressure, \( p_{star} \) [Pa]
        real(real64) :: p_star
        !> Specific Gas Constant, \( R \) [J/(kg K)]
        real(real64) :: R
        !> Initialization flag
        logical :: is_initialized = .false.
    contains
        procedure(abst_initialize_iapws_gibbs), pass(self), public, deferred :: initialize
        procedure(abst_calc_gibbs_iapws), pass(self), public, deferred :: calc_gamma
        procedure, pass(self), public :: calc_properties => calc_properties_gibbs
        procedure, pass(self), public :: calc_nu => calc_nu_gibbs
        procedure, pass(self), public :: calc_rho => calc_rho_gibbs
        procedure, pass(self), public :: calc_drho_dT => calc_drho_dT_gibbs
        procedure, pass(self), public :: calc_drho_dP => calc_drho_dP_gibbs
        procedure, pass(self), public :: calc_u => calc_u_gibbs
        procedure, pass(self), public :: calc_s => calc_s_gibbs
        procedure, pass(self), public :: calc_h => calc_h_gibbs
        procedure, pass(self), public :: calc_cp => calc_cp_gibbs
        procedure, pass(self), public :: calc_cv => calc_cv_gibbs
        procedure, pass(self), public :: calc_w => calc_w_gibbs
        procedure, pass(self), public :: calc_alpha => calc_alpha_gibbs
        procedure, pass(self), public :: calc_beta => calc_beta_gibbs
        procedure, pass(self), public :: calc_kappa_T => calc_kappa_T_gibbs
        procedure, pass(self), public :: calc_kappa_s => calc_kappa_s_gibbs
    end type abst_iapws_gibbs

    !> Data structure to hold IAPWS properties
    type :: type_iapws_helmholtz_property
        !> Total adimensional helmholtz energy [-]
        real(real64) :: f
        !> Derivative of total adimensional helmholtz energy with respect to delta,∂phi/∂δ|τ [-]
        real(real64) :: f_d
        !> Derivative of total adimensional helmholtz energy with respect to tau,∂phi/∂τ|δ [-]
        real(real64) :: f_t
        !> Second derivative of total adimensional helmholtz energy with respect to delta,∂²phi/∂δ²|τ [-]
        real(real64) :: f_dd
        !> Second derivative of total adimensional helmholtz energy with respect to tau,∂²phi/∂τ²|δ [-]
        real(real64) :: f_tt
        !> Mixed second derivative of total adimensional helmholtz energy,∂²phi/∂δ∂τ [-]
        real(real64) :: f_dt
        !> Ideal adimensional helmholtz energy [-]
        real(real64) :: f0
        !> Ideal adimensional helmholtz energy derivative with respect to delta,∂phi0/∂δ|τ [-]
        real(real64) :: f0_d
        !> Ideal adimensional helmholtz energy derivative with respect to tau,∂phi0/∂τ|δ [-]
        real(real64) :: f0_t
        !> Ideal adimensional helmholtz energy second derivative with respect to delta,∂²phi0/∂δ²|τ [-]
        real(real64) :: f0_dd
        !> Ideal adimensional helmholtz energy second derivative with respect to tau,∂²phi0/∂τ²|δ [-]
        real(real64) :: f0_tt
        !> Ideal adimensional helmholtz energy mixed second derivative,∂²phi0/∂δ∂τ [-]
        real(real64) :: f0_dt
        !> Ideal adimensional helmholtz energy [-]
        real(real64) :: fr
        !> Ideal adimensional helmholtz energy derivative with respect to delta,∂phir/∂δ|τ [-]
        real(real64) :: fr_d
        !> Ideal adimensional helmholtz energy derivative with respect to tau,∂phir/∂τ|δ [-]
        real(real64) :: fr_t
        !> Ideal adimensional helmholtz energy second derivative with respect to delta,∂²phir/∂δ²|τ [-]
        real(real64) :: fr_dd
        !> Ideal adimensional helmholtz energy second derivative with respect to tau,∂²phir/∂τ²|δ [-]
        real(real64) :: fr_tt
        !> Ideal adimensional helmholtz energy mixed second derivative,∂²phir/∂δ∂τ [-]
        real(real64) :: fr_dt
    contains
        procedure, pass(self), public :: reset => reset_iapws_helmholtz_property
    end type type_iapws_helmholtz_property

    !> Data structure to hold IAPWS Gibbs coefficients
    type :: type_iapws_gibbs_coefficient
        !> Total adimensional gibbs energy [-]
        real(real64) :: g
        !> Derivative of total adimensional gibbs energy with respect to pi,∂γ/∂π|τ [-]
        real(real64) :: g_p
        !> Derivative of total adimensional gibbs energy with respect to tau,∂γ/∂τ|π [-]
        real(real64) :: g_t
        !> Second derivative of total adimensional gibbs energy with respect to pi,∂²γ/∂π²|τ [-]
        real(real64) :: g_pp
        !> Second derivative of total adimensional gibbs energy with respect to tau,∂²γ/∂τ²|π [-]
        real(real64) :: g_tt
        !> Mixed second derivative of total adimensional gibbs energy,∂²γ/∂π∂τ [-]
        real(real64) :: g_pt
        !> Ideal adimensional gibbs energy [-]
        real(real64) :: g0
        !> Ideal adimensional gibbs energy derivative with respect to pi,∂γ0/∂π|τ [-]
        real(real64) :: g0_p
        !> Ideal adimensional gibbs energy derivative with respect to tau,∂γ0/∂τ|π [-]
        real(real64) :: g0_t
        !> Ideal adimensional gibbs energy second derivative with respect to pi,∂²γ0/∂π²|τ [-]
        real(real64) :: g0_pp
        !> Ideal adimensional gibbs energy second derivative with respect to tau,∂²γ0/∂τ²|π [-]
        real(real64) :: g0_tt
        !> Ideal adimensional gibbs energy mixed second derivative,∂²γ0/∂π∂τ [-]
        real(real64) :: g0_pt
        !> Residual adimensional gibbs energy [-]
        real(real64) :: gr
        !> Residual adimensional gibbs energy derivative with respect to pi,∂γr/∂π|τ [-]
        real(real64) :: gr_p
        !> Residual adimensional gibbs energy derivative with respect to tau,∂γr/∂τ|π [-]
        real(real64) :: gr_t
        !> Residual adimensional gibbs energy second derivative with respect to pi,∂²γr/∂π²|τ [-]
        real(real64) :: gr_pp
        !> Residual adimensional gibbs energy second derivative with respect to tau,∂²γr/∂τ²|π [-]
        real(real64) :: gr_tt
        !> Residual adimensional gibbs energy mixed second derivative,∂²γr/∂π∂τ [-]
        real(real64) :: gr_pt
    contains
        procedure, pass(self), public :: reset => reset_iapws_gibbs_coefficient
    end type type_iapws_gibbs_coefficient

    interface
        module pure elemental subroutine reset_iapws_helmholtz_property(self)
            implicit none
            !> IAPWS Helmholtz properties object
            class(type_iapws_helmholtz_property), intent(inout) :: self
        end subroutine reset_iapws_helmholtz_property

        module pure elemental subroutine reset_iapws_gibbs_coefficient(self)
            implicit none
            !> IAPWS Gibbs properties object
            class(type_iapws_gibbs_coefficient), intent(inout) :: self
        end subroutine reset_iapws_gibbs_coefficient
    end interface

    !> IAPWS properties data structure
    type :: type_iapws_property
        !> Region ID
        integer(int32) :: region_id = -1
        !> Specific Volume, \( \nu \) [m^3/kg]
        real(real64) :: nu = 0.0d0
        !> Density, \( \rho \) [kg/m^3]
        real(real64) :: rho = 0.0d0
        !> Temperature derivative of density, \( \left(\frac{\partial \rho}{\partial T}\right)_P \) [kg/(m^3 K)]
        real(real64) :: drho_dT = 0.0d0
        !> Pressure derivative of density, \( \left(\frac{\partial \rho}{\partial P}\right)_T \) [kg/(m^3 Pa)]
        real(real64) :: drho_dP = 0.0d0
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
            import :: abst_iapws_helmholtz, type_iapws_helmholtz_property, real64
            implicit none
            !> IAPWS Helmholtz model instance
            class(abst_iapws_helmholtz), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced density ρ/rho_c, [-]
            real(real64), intent(in) :: delta
            !> IAPWS-95 ideal helmholtz properties
            type(type_iapws_helmholtz_property), intent(inout) :: property
        end subroutine abst_calc_phi_iapws
    end interface

    abstract interface
        pure elemental subroutine abst_initialize_iapws_gibbs(self)
            import :: abst_iapws_gibbs
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(inout) :: self
        end subroutine abst_initialize_iapws_gibbs

        !> Calculate IAPWS Gibbs coefficients with given T and P
        pure elemental subroutine abst_calc_gibbs_iapws(self, T_in, P_in, coef)
            import :: abst_iapws_gibbs, type_iapws_gibbs_coefficient, real64
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: P_in
            !> IAPWS Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        end subroutine abst_calc_gibbs_iapws
    end interface

    interface
        module pure elemental subroutine calc_properties_helmholtz(self, T_in, rho_in, property)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            type(type_iapws_property), intent(inout) :: property

        end subroutine calc_properties_helmholtz

        module pure elemental subroutine calc_rho_helmholtz(self, T_in, P_in, rho)
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: P_in
            real(real64), intent(inout) :: rho

        end subroutine calc_rho_helmholtz

        module pure elemental subroutine calc_p_helmholtz(self, T_in, rho_in, p, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: p
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_p_helmholtz

        module pure elemental subroutine calc_p_rho_helmholtz(self, T_in, rho_in, p_rho, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: p_rho
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_p_rho_helmholtz

        module pure elemental subroutine calc_p_T_helmholtz(self, T_in, rho_in, p_T, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: p_T
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_p_T_helmholtz

        module pure elemental subroutine calc_u_helmholtz(self, T_in, rho_in, u, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: u
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_u_helmholtz

        module pure elemental subroutine calc_s_helmholtz(self, T_in, rho_in, s, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: s
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_s_helmholtz

        module pure elemental subroutine calc_h_helmholtz(self, T_in, rho_in, h, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: h
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_h_helmholtz

        module pure elemental subroutine calc_cp_helmholtz(self, T_in, rho_in, cp, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: cp
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_cp_helmholtz

        module pure elemental subroutine calc_cv_helmholtz(self, T_in, rho_in, cv, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: cv
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_cv_helmholtz

        module pure elemental subroutine calc_w_helmholtz(self, T_in, rho_in, w, prop_in)
            implicit none
            class(abst_iapws_helmholtz), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: rho_in
            real(real64), intent(inout) :: w
            type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        end subroutine calc_w_helmholtz
    end interface

    interface
        module pure elemental subroutine calc_properties_gibbs(self, T_in, p_in, property)
            implicit none
            class(abst_iapws_gibbs), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            type(type_iapws_property), intent(inout) :: property

        end subroutine calc_properties_gibbs

        !> Calculate specific volume using Gibbs formulation
        module pure elemental subroutine calc_nu_gibbs(self, T_in, p_in, nu, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific volume, ν [m^3/kg].
            real(real64), intent(inout) :: nu
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_nu_gibbs

        !> Calculate density using Gibbs formulation
        module pure elemental subroutine calc_rho_gibbs(self, T_in, p_in, rho, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Density, ρ [kg/m^3].
            real(real64), intent(inout) :: rho
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_rho_gibbs

        !> 密度の温度微分 (drho/dT)_P
        module pure elemental subroutine calc_drho_dT_gibbs(self, T_in, P_in, drho_dT, prop_in)
            implicit none
            class(abst_iapws_gibbs), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: P_in
            real(real64), intent(inout) :: drho_dT
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_drho_dT_gibbs

        !> 密度の圧力微分 (drho/dP)_T
        module pure elemental subroutine calc_drho_dP_gibbs(self, T_in, P_in, drho_dP, prop_in)
            implicit none
            class(abst_iapws_gibbs), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: P_in
            real(real64), intent(inout) :: drho_dP
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_drho_dP_gibbs

        !> Calculate specific internal energy using Gibbs formulation
        module pure elemental subroutine calc_u_gibbs(self, T_in, p_in, u, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific internal energy, u [J/kg].
            real(real64), intent(inout) :: u
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_u_gibbs

        !> Calculate specific entropy using Gibbs formulation
        module pure elemental subroutine calc_s_gibbs(self, T_in, p_in, s, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific entropy, s [J/(kg·K)].
            real(real64), intent(inout) :: s
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_s_gibbs

        !> Calculate specific enthalpy using Gibbs formulation
        module pure elemental subroutine calc_h_gibbs(self, T_in, p_in, h, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific enthalpy, h [J/kg].
            real(real64), intent(inout) :: h
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_h_gibbs

        !> Calculate specific heat capacity at constant pressure using Gibbs formulation
        module pure elemental subroutine calc_cp_gibbs(self, T_in, p_in, cp, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific heat capacity at constant pressure, cp [J/(kg·K)].
            real(real64), intent(inout) :: cp
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_cp_gibbs

        !> Calculate specific heat capacity at constant volume using Gibbs formulation
        module pure elemental subroutine calc_cv_gibbs(self, T_in, p_in, cv, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Specific heat capacity at constant volume, cv [J/(kg·K)].
            real(real64), intent(inout) :: cv
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_cv_gibbs

        !> Calculate speed of sound using Gibbs formulation
        module pure elemental subroutine calc_w_gibbs(self, T_in, p_in, w, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: p_in
            !> Speed of sound, w [m/s].
            real(real64), intent(inout) :: w
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_w_gibbs

        !> Calculate cubic expansion coefficient using Gibbs formulation
        module pure elemental subroutine calc_alpha_gibbs(self, T_in, P_in, alpha, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: P_in
            !> Cubic expansion coefficient [1/K].
            real(real64), intent(inout) :: alpha
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in
        end subroutine calc_alpha_gibbs

        !> Calculate pressure coefficient using Gibbs formulation
        module pure elemental subroutine calc_beta_gibbs(self, T_in, P_in, beta, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: P_in
            !> Pressure coefficient [Pa/K].
            real(real64), intent(inout) :: beta
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_beta_gibbs

        !> Calculate isothermal compressibility using Gibbs formulation
        module pure elemental subroutine calc_kappa_T_gibbs(self, T_in, P_in, kappa_T, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: P_in
            !> Isothermal compressibility [1/Pa].
            real(real64), intent(inout) :: kappa_T
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_kappa_T_gibbs

        !> Calculate adiabatic compressibility using Gibbs formulation
        module pure elemental subroutine calc_kappa_s_gibbs(self, T_in, P_in, kappa_s, prop_in)
            implicit none
            !> IAPWS Gibbs model instance
            class(abst_iapws_gibbs), intent(in) :: self
            !> Temperature, T [K].
            real(real64), intent(in) :: T_in
            !> Pressure, P [Pa].
            real(real64), intent(in) :: P_in
            !> Adiabatic compressibility [1/Pa].
            real(real64), intent(inout) :: kappa_s
            !> Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        end subroutine calc_kappa_s_gibbs
    end interface

end module module_iapws
