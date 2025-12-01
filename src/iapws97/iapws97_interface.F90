module module_iapws97
    use, intrinsic :: iso_fortran_env
    use :: module_kahan
    use :: iapws
    implicit none

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
        module pure elemental subroutine calc_gamma_iapws97_region1(self, tau, pi, property)
            implicit none
            class(type_iapws97_region1), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced pressure p/p_c, [-]
            real(real64), intent(in) :: pi
            !> IAPWS Gibbs properties
            type(type_iapws_gamma_property), intent(inout) :: property
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

        !> Calculate the dimensionless Gibbs free energy \(\gamma\) for Region 2.
        module pure elemental subroutine calc_gamma_iapws97_region2(self, tau, pi, property)
            implicit none
            !> IAPWS-97 Region 2 object.
            class(type_iapws97_region2), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced pressure p/p_c, [-]
            real(real64), intent(in) :: pi
            !> IAPWS Gibbs properties
            type(type_iapws_gamma_property), intent(inout) :: property
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
            type(type_iapws_phi_property), intent(inout) :: property
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
        module pure elemental subroutine calc_gamma_iapws97_region5(self, tau, pi, property)
            implicit none
            !> IAPWS-97 Region 5 object.
            class(type_iapws97_region5), intent(in) :: self
            !> Inverse reduced temperature Tc/T, [-]
            real(real64), intent(in) :: tau
            !> Reduced pressure p/p_c, [-]
            real(real64), intent(in) :: pi
            !> IAPWS Gibbs properties
            type(type_iapws_gamma_property), intent(inout) :: property
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
        procedure, pass(self), public :: calc_u => calc_u_iapws97
        procedure, pass(self), public :: calc_h => calc_h_iapws97
        procedure, pass(self), public :: calc_s => calc_s_iapws97
        procedure, pass(self), public :: calc_cp => calc_cp_iapws97
        procedure, pass(self), public :: calc_cv => calc_cv_iapws97
        procedure, pass(self), public :: calc_w => calc_w_iapws97
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

        module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: nu
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_nu_iapws97

        module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: rho
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_rho_iapws97

        module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: u
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_u_iapws97

        module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: h
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_h_iapws97

        module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: s
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_s_iapws97

        module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: cp
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_cp_iapws97

        module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: cv
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_cv_iapws97

        module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w, prop_in)
            implicit none
            class(type_iapws97), intent(in) :: self
            real(real64), intent(in) :: T_in
            real(real64), intent(in) :: p_in
            real(real64), intent(inout) :: w
            type(type_iapws_property), intent(inout), optional :: prop_in

        end subroutine calc_w_iapws97
    end interface

end module module_iapws97
