module module_iapws06_iceIh
    use, intrinsic :: iso_fortran_env
    use :: module_iapws
    implicit none

    type, extends(abst_iapws_gibbs) :: type_iapws06
        private
        real(real64) :: pi_0
    contains
        procedure, pass(self), public :: initialize => initialize_type_iapws06
        procedure, pass(self), public :: calc_gamma => calc_gamma_iapws06
        procedure, pass(self), public :: calc_dcp_dT => calc_dcp_dT_iapws06
        procedure, pass(self), public :: calc_dcp_dP => calc_dcp_dP_iapws06
        procedure, pass(self), private :: calc_gamma_Ih => calc_gamma_iapws06_Ih
        procedure, pass(self), private :: calc_gamma_t_Ih => calc_gamma_t_iapws06_Ih
        procedure, pass(self), private :: calc_gamma_p_Ih => calc_gamma_p_iapws06_Ih
        procedure, pass(self), private :: calc_gamma_tt_Ih => calc_gamma_tt_iapws06_Ih
        procedure, pass(self), private :: calc_gamma_tp_Ih => calc_gamma_tp_iapws06_Ih
        procedure, pass(self), private :: calc_gamma_pp_Ih => calc_gamma_pp_iapws06_Ih
    end type type_iapws06

    interface
        !> Isobaric heat-capacity temperature derivative, [J/(kg K^2)].
        !> Returns quiet NaN outside 0 < T <= 273.16 K, 0 <= P <= 210 MPa.
        module pure elemental subroutine calc_dcp_dT_iapws06(self, T_in, P_in, dcp_dT)
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: T_in, P_in
            real(real64), intent(out) :: dcp_dT
        end subroutine calc_dcp_dT_iapws06

        !> Isothermal heat-capacity pressure derivative, [J/(kg K Pa)].
        !> Returns quiet NaN outside 0 < T <= 273.16 K, 0 <= P <= 210 MPa.
        module pure elemental subroutine calc_dcp_dP_iapws06(self, T_in, P_in, dcp_dP)
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: T_in, P_in
            real(real64), intent(out) :: dcp_dP
        end subroutine calc_dcp_dP_iapws06

        !> Initialize the IAPWS-06 object.
        module pure elemental subroutine initialize_type_iapws06(self)
            implicit none
            !> Initialize the IAPWS-06 object.
            class(type_iapws06), intent(inout) :: self
        end subroutine initialize_type_iapws06

        !> Calculate the dimensionless Gibbs free energy \(\gamma\).
        module pure elemental subroutine calc_gamma_iapws06(self, T_in, P_in, coef)
            implicit none
            !> Initialize the IAPWS-06 object.
            class(type_iapws06), intent(in) :: self
            !> Input temperature, [K]
            real(real64), intent(in) :: T_in
            !> Input pressure, [Pa]
            real(real64), intent(in) :: P_in
            !> IAPWS Gibbs coefficients
            type(type_iapws_gibbs_coefficient), intent(inout) :: coef

        end subroutine calc_gamma_iapws06

        module pure elemental function calc_gamma_iapws06_Ih(self, tau, pi) result(gamma)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma

        end function calc_gamma_iapws06_Ih

        module pure elemental function calc_gamma_t_iapws06_Ih(self, tau, pi) result(gamma_t)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma_t

        end function calc_gamma_t_iapws06_Ih

        module pure elemental function calc_gamma_p_iapws06_Ih(self, tau, pi) result(gamma_p)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma_p

        end function calc_gamma_p_iapws06_Ih

        module pure elemental function calc_gamma_tt_iapws06_Ih(self, tau, pi) result(gamma_tt)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma_tt

        end function calc_gamma_tt_iapws06_Ih

        module pure elemental function calc_gamma_tp_iapws06_Ih(self, tau, pi) result(gamma_tp)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma_tp

        end function calc_gamma_tp_iapws06_Ih

        module pure elemental function calc_gamma_pp_iapws06_Ih(self, tau, pi) result(gamma_pp)
            implicit none
            class(type_iapws06), intent(in) :: self
            real(real64), intent(in) :: tau
            real(real64), intent(in) :: pi
            real(real64) :: gamma_pp

        end function calc_gamma_pp_iapws06_Ih
    end interface

end module module_iapws06_iceIh
