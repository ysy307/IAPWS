submodule(module_iapws06_iceIh) iapws06_iceIh_base
    implicit none

    !---------------------------------------------------------------------------
    ! Constants for IAPWS-06 Ice Ih
    !---------------------------------------------------------------------------
    ! Triple point pressure used for scaling [Pa]
    real(real64), parameter :: p_starIh = 611.657d0

    real(real64), parameter :: s0_Ih = -0.332733756492168d4
    integer(int32), parameter :: gamma0_terms = 5
    real(real64), parameter :: gamma0_Ih(gamma0_terms) = [ &
                               -0.632020233335886d6, &
                               0.655022213658955d0, &
                               -0.189369929326131d-7, &
                               0.339746123271053d-14, &
                               -0.556464869058991d-21]
    complex(real64), parameter :: t1_Ih = (0.368017112855051d-1, 0.510878114959572d-1)
    complex(real64), parameter :: r1_Ih = (0.447050716285388d2, 0.656876847463481d2)
    complex(real64), parameter :: t2_Ih = (0.337315741065416d0, 0.335449415919309d0)
    integer(int32), parameter :: r2_terms = 3
    complex(real64), parameter :: r2_Ih_r2(r2_terms) = [ &
                                  (-0.725974574329220d2, -0.781008427112870d2), &
                                  (-0.557107698030123d-4, 0.464578634580806d-4), &
                                  (0.234801409215913d-10, -0.285651142904972d-10)]

contains

    !> Initialize the IAPWS-06 object.
    module pure elemental subroutine initialize_type_iapws06(self)
        implicit none
        class(type_iapws06), intent(inout) :: self

        self%T_star = 273.16d0 ! [K] Triple point temperature
        self%p_star = p_starIh ! [Pa] Triple point pressure

        ! Reference dimensionless pressure pi_0 corresponding to P = 101325 Pa
        self%pi_0 = 101325.0d0 / self%p_star

        self%R = 461.526d0 ! [J/kgK]
    end subroutine initialize_type_iapws06

    !> Calculate the Gibbs free energy properties for Ice Ih.
    module pure elemental subroutine calc_gamma_iapws06(self, T_in, P_in, coef)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: T_in ! [K]
        real(real64), intent(in) :: P_in ! [Pa]
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef

        real(real64) :: tau, pi

        ! IAPWS-06 uses tau = T / Tt (Direct ratio, unlike IAPWS-97)
        tau = T_in / self%T_star
        ! pi = P / Pt
        pi = P_in / self%p_star

        call coef%reset()

        ! Calculate properties and store in coef
        ! These functions return Dimensional values (J/kg, m3/kg, etc.)
        coef%g = calc_gamma_iapws06_Ih(self, tau, pi)
        coef%g_p = calc_gamma_p_iapws06_Ih(self, tau, pi)
        coef%g_t = calc_gamma_t_iapws06_Ih(self, tau, pi)
        coef%g_pp = calc_gamma_pp_iapws06_Ih(self, tau, pi)
        coef%g_tt = calc_gamma_tt_iapws06_Ih(self, tau, pi)
        coef%g_pt = calc_gamma_tp_iapws06_Ih(self, tau, pi)

    end subroutine calc_gamma_iapws06

    !---------------------------------------------------------------------------
    ! Gibbs Energy: g(T,p) [J/kg]
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_iapws06_Ih(self, tau, pi) result(gamma)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: tau, pi
        real(real64) :: gamma

        complex(real64) :: r2
        complex(real64) :: term1, term2

        ! Note: Passing self%pi_0 to helpers
        gamma = calc_gamma0_iapws06_Ih(pi, self%pi_0) - s0_Ih * self%T_star * tau
        r2 = calc_r2_iapws06_Ih(pi, self%pi_0)

        ! k=1
        term1 = (t1_Ih - tau) * log(t1_Ih - tau) + (t1_Ih + tau) * log(t1_Ih + tau) &
                - 2.0d0 * t1_Ih * log(t1_Ih) - (tau**2.0d0) / t1_Ih

        ! k=2
        term2 = (t2_Ih - tau) * log(t2_Ih - tau) + (t2_Ih + tau) * log(t2_Ih + tau) &
                - 2.0d0 * t2_Ih * log(t2_Ih) - (tau**2.0d0) / t2_Ih

        gamma = gamma + self%T_star * dble(r1_Ih * term1 + r2 * term2)
    end function calc_gamma_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Derivative w.r.t Temperature: g_T [J/(kg K)] = -Entropy
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_t_iapws06_Ih(self, tau, pi) result(gamma_t)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: pi, tau
        real(real64) :: gamma_t

        complex(real64) :: r2
        complex(real64) :: term1, term2

        r2 = calc_r2_iapws06_Ih(pi, self%pi_0)

        ! k=1
        term1 = -log(t1_Ih - tau) + log(t1_Ih + tau) - 2.0d0 * tau / t1_Ih
        ! k=2
        term2 = -log(t2_Ih - tau) + log(t2_Ih + tau) - 2.0d0 * tau / t2_Ih

        gamma_t = -s0_Ih + dble(r1_Ih * term1 + r2 * term2)
    end function calc_gamma_t_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Derivative w.r.t Pressure: g_p [m^3/kg] = Specific Volume
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_p_iapws06_Ih(self, tau, pi) result(gamma_p)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: tau, pi
        real(real64) :: gamma_p

        complex(real64) :: r2_p
        complex(real64) :: term2

        r2_p = calc_r2_p_iapws06_Ih(pi, self%pi_0)

        ! k=2 only
        term2 = (t2_Ih - tau) * log(t2_Ih - tau) + (t2_Ih + tau) * log(t2_Ih + tau) &
                - 2.0d0 * t2_Ih * log(t2_Ih) - (tau**2.0d0) / t2_Ih

        gamma_p = calc_gamma0_p_iapws06_Ih(pi, self%pi_0) + self%T_star * dble(r2_p * term2)
    end function calc_gamma_p_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Second Derivative w.r.t Temperature: g_TT [J/(kg K^2)]
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_tt_iapws06_Ih(self, tau, pi) result(gamma_tt)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: tau, pi
        real(real64) :: gamma_tt

        complex(real64) :: r2
        complex(real64) :: term1, term2

        r2 = calc_r2_iapws06_Ih(pi, self%pi_0)

        ! k=1
        term1 = 1.0d0 / (t1_Ih - tau) + 1.0d0 / (t1_Ih + tau) - 2.0d0 / t1_Ih
        ! k=2
        term2 = 1.0d0 / (t2_Ih - tau) + 1.0d0 / (t2_Ih + tau) - 2.0d0 / t2_Ih

        gamma_tt = (1.0d0 / self%T_star) * dble(r1_Ih * term1 + r2 * term2)
    end function calc_gamma_tt_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Mixed Derivative: g_Tp [m^3/(kg K)]
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_tp_iapws06_Ih(self, tau, pi) result(gamma_tp)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: tau, pi
        real(real64) :: gamma_tp

        complex(real64) :: r2_p
        complex(real64) :: term2

        r2_p = calc_r2_p_iapws06_Ih(pi, self%pi_0)

        ! k=2 only
        term2 = -log(t2_Ih - tau) + log(t2_Ih + tau) - 2.0d0 * tau / t2_Ih

        gamma_tp = dble(r2_p * term2)
    end function calc_gamma_tp_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Second Derivative w.r.t Pressure: g_pp [m^3/(kg Pa)]
    !---------------------------------------------------------------------------
    module pure elemental function calc_gamma_pp_iapws06_Ih(self, tau, pi) result(gamma_pp)
        implicit none
        class(type_iapws06), intent(in) :: self
        real(real64), intent(in) :: tau, pi
        real(real64) :: gamma_pp

        complex(real64) :: r2_pp
        complex(real64) :: term2

        r2_pp = calc_r2_pp_iapws06_Ih(pi)

        ! k=2 only
        term2 = (t2_Ih - tau) * log(t2_Ih - tau) + (t2_Ih + tau) * log(t2_Ih + tau) &
                - 2.0d0 * t2_Ih * log(t2_Ih) - (tau**2.0d0) / t2_Ih

        gamma_pp = calc_gamma0_pp_iapws06_Ih(pi, self%pi_0) + self%T_star * dble(r2_pp * term2)
    end function calc_gamma_pp_iapws06_Ih

    !---------------------------------------------------------------------------
    ! Helper Functions
    !---------------------------------------------------------------------------

    pure elemental function calc_gamma0_iapws06_Ih(pi, pi_0) result(gamma0)
        implicit none
        real(real64), intent(in) :: pi, pi_0
        real(real64) :: gamma0
        real(real64) :: delta_pi
        integer(int32) :: i

        delta_pi = pi - pi_0
        gamma0 = gamma0_Ih(gamma0_terms)
        do i = gamma0_terms - 1, 1, -1
            gamma0 = gamma0_Ih(i) + delta_pi * gamma0
        end do
    end function calc_gamma0_iapws06_Ih

    pure elemental function calc_gamma0_p_iapws06_Ih(pi, pi_0) result(val)
        implicit none
        real(real64), intent(in) :: pi, pi_0
        real(real64) :: val, delta_pi
        integer(int32) :: i

        delta_pi = pi - pi_0
        val = gamma0_Ih(5) * 4.0d0
        do i = 4, 2, -1
            val = gamma0_Ih(i) * dble(i - 1) + delta_pi * val
        end do
        ! Use module parameter p_starIh (611.657)
        val = val / p_starIh
    end function calc_gamma0_p_iapws06_Ih

    pure elemental function calc_gamma0_pp_iapws06_Ih(pi, pi_0) result(val)
        implicit none
        real(real64), intent(in) :: pi, pi_0
        real(real64) :: val, delta_pi
        integer(int32) :: i

        delta_pi = pi - pi_0
        val = gamma0_Ih(5) * 12.0d0
        do i = 4, 3, -1
            val = gamma0_Ih(i) * dble((i - 1) * (i - 2)) + delta_pi * val
        end do
        ! Use module parameter p_starIh
        val = val / (p_starIh**2.0d0)
    end function calc_gamma0_pp_iapws06_Ih

    pure elemental function calc_r2_iapws06_Ih(pi, pi_0) result(val)
        implicit none
        real(real64), intent(in) :: pi, pi_0
        complex(real64) :: val
        real(real64) :: delta_pi
        integer(int32) :: i

        delta_pi = pi - pi_0
        val = r2_Ih_r2(r2_terms)
        do i = r2_terms - 1, 1, -1
            val = r2_Ih_r2(i) + delta_pi * val
        end do
    end function calc_r2_iapws06_Ih

    pure elemental function calc_r2_p_iapws06_Ih(pi, pi_0) result(val)
        implicit none
        real(real64), intent(in) :: pi, pi_0
        complex(real64) :: val
        real(real64) :: delta_pi
        integer(int32) :: i

        delta_pi = pi - pi_0
        val = r2_Ih_r2(3) * 2.0d0
        do i = 2, 2, -1
            val = r2_Ih_r2(i) * dble(i - 1) + delta_pi * val
        end do
        val = val / p_starIh
    end function calc_r2_p_iapws06_Ih

    pure elemental function calc_r2_pp_iapws06_Ih(pi) result(val)
        implicit none
        real(real64), intent(in) :: pi
        complex(real64) :: val
        val = r2_Ih_r2(3) * 2.0d0 / (p_starIh**2.0d0)
    end function calc_r2_pp_iapws06_Ih

end submodule iapws06_iceIh_base
