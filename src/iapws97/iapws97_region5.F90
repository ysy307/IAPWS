submodule(module_iapws97) iapws97_region5
    implicit none
    !------------------------------------------------------------------------------------------
    ! Region5: High temperature region (IAPWS-IF97)
    !------------------------------------------------------------------------------------------
    integer(int32), parameter :: N05_terms = 6
    integer(int32), parameter :: J0_r5(N05_terms) = [0, 1, -3, -2, -1, 2]
    real(real64), parameter :: n0_r5(N05_terms) = [ &
                               -0.13179983674201d2, 0.68540841634434d1, -0.24805148933466d-1, &
                               0.36901534980333d0, -0.31161318213925d1, -0.32961626538917d0]

    integer(int32), parameter :: Nr5_terms = 6
    integer(int32), parameter :: Ir_r5(Nr5_terms) = [1, 1, 1, 2, 2, 3]
    integer(int32), parameter :: Jr_r5(Nr5_terms) = [1, 2, 3, 3, 9, 7]
    real(real64), parameter :: nr_r5(Nr5_terms) = [ &
                               0.15736404855259d-2, 0.90153761673944d-3, -0.50270077677648d-2, &
                               0.22440037409485d-5, -0.41163275453471d-5, 0.37919454822955d-7]
contains

    !> Initialize the IAPWS-97 Region 5 object.
    module pure elemental subroutine initialize_type_iapws97_region5(self)
        implicit none
        !> IAPWS-97 Region 5 object.
        class(type_iapws97_region5), intent(inout) :: self

        self%T_star = 1000.0d0 ! [K]
        self%p_star = 1.0d6 ! [Pa]
        self%R = 461.526d0 ! [J/(kg·K)]

        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_region5

    !> Calculate the Reduced Gibbs free energy \(\gamma\) for Region 5.
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

        call property%reset()

        ! Calculate gamma and its derivatives
        property%gamma0 = calc_gamma0_region5(tau, pi)
        property%gamma0_p = calc_gamma0_p_region5(tau, pi)
        property%gamma0_t = calc_gamma0_t_region5(tau, pi)
        property%gamma0_pp = calc_gamma0_pp_region5(tau, pi)
        property%gamma0_tt = calc_gamma0_tt_region5(tau, pi)
        property%gamma0_pt = calc_gamma0_pt_region5(tau, pi)
        property%gammar = calc_gammar_region5(tau, pi)
        property%gammar_p = calc_gammar_p_region5(tau, pi)
        property%gammar_t = calc_gammar_t_region5(tau, pi)
        property%gammar_pp = calc_gammar_pp_region5(tau, pi)
        property%gammar_tt = calc_gammar_tt_region5(tau, pi)
        property%gammar_pt = calc_gammar_pt_region5(tau, pi)

        property%gamma = property%gamma0 + property%gammar
        property%gamma_p = property%gamma0_p + property%gammar_p
        property%gamma_t = property%gamma0_t + property%gammar_t
        property%gamma_pp = property%gamma0_pp + property%gammar_pp
        property%gamma_tt = property%gamma0_tt + property%gammar_tt
        property%gamma_pt = property%gamma0_pt + property%gammar_pt
    end subroutine calc_gamma_iapws97_region5

    pure elemental function calc_gamma0_region5(tau, pi) result(gamma0)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Ideal gas part of reduced Gibbs free energy, \(\gamma_0\)
        real(real64) :: gamma0

        integer(int32) :: i
        real(real64) :: c_gamma

        gamma0 = log(pi)
        c_gamma = 0.0d0

        do i = 1, N05_terms
            call kahan_add(gamma0, c_gamma, &
                           n0_r5(i) * (tau**J0_r5(i)))
        end do
    end function calc_gamma0_region5

    pure elemental function calc_gamma0_p_region5(tau, pi) result(gamma0_p)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> First derivative \(\gamma^0_{\pi}\)
        real(real64) :: gamma0_p

        gamma0_p = 1.0d0 / pi
    end function calc_gamma0_p_region5

    pure elemental function calc_gamma0_pp_region5(tau, pi) result(gamma0_pp)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Second derivative \(\gamma^0_{\pi \pi}\)
        real(real64) :: gamma0_pp

        ! Table 40: Second derivative w.r.t pi is -1/pi^2
        gamma0_pp = -1.0d0 / (pi * pi)
    end function calc_gamma0_pp_region5

    pure elemental function calc_gamma0_t_region5(tau, pi) result(gamma0_t)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> First derivative \(\gamma^0_{t}\)
        real(real64) :: gamma0_t

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma0_t

        gamma0_t = 0.0d0
        c_gamma0_t = 0.0d0
        do i = 1, N05_terms
            J_val = real(J0_r5(i), real64)
            call kahan_add(gamma0_t, c_gamma0_t, &
                           n0_r5(i) * J_val * (tau**(J0_r5(i) - 1)))
        end do
    end function calc_gamma0_t_region5

    pure elemental function calc_gamma0_tt_region5(tau, pi) result(gamma0_tt)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Second derivative \(\gamma^0_{tt}\)
        real(real64) :: gamma0_tt

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma0_tt

        gamma0_tt = 0.0d0
        c_gamma0_tt = 0.0d0

        do i = 1, N05_terms
            J_val = real(J0_r5(i), real64)
            call kahan_add(gamma0_tt, c_gamma0_tt, &
                           n0_r5(i) * J_val * (J_val - 1.0d0) * (tau**(J0_r5(i) - 2)))
        end do
    end function calc_gamma0_tt_region5

    pure elemental function calc_gamma0_pt_region5(tau, pi) result(gamma0_pt)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Mixed derivative \(\gamma^0_{p t}\)
        real(real64) :: gamma0_pt

        gamma0_pt = 0.0d0
    end function calc_gamma0_pt_region5

    pure elemental function calc_gammar_region5(tau, pi) result(gammar)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Helmholtz free energy residual part \(\gamma^r\)
        real(real64) :: gammar

        integer(int32) :: i
        real(real64) :: c_gammar

        gammar = 0.0d0
        c_gammar = 0.0d0

        do i = 1, Nr5_terms
            call kahan_add(gammar, c_gammar, &
                           nr_r5(i) * (pi**Ir_r5(i)) * (tau**Jr_r5(i)))
        end do
    end function calc_gammar_region5

    pure elemental function calc_gammar_p_region5(tau, pi) result(gammar_p)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> First derivative \(\gamma^r_{\pi}\)
        real(real64) :: gammar_p

        integer(int32) :: i
        integer(int32) :: I_val
        real(real64) :: c_gammar_p

        gammar_p = 0.0d0
        c_gammar_p = 0.0d0
        do i = 1, Nr5_terms
            I_val = real(Ir_r5(i), real64)
            call kahan_add(gammar_p, c_gammar_p, &
                           nr_r5(i) * I_val * (pi**(Ir_r5(i) - 1)) * (tau**Jr_r5(i)))
        end do
    end function calc_gammar_p_region5

    pure elemental function calc_gammar_pp_region5(tau, pi) result(gammar_pp)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Second derivative \(\gamma^r_{\pi \pi}\)
        real(real64) :: gammar_pp

        integer(int32) :: i
        integer(int32) :: I_val
        real(real64) :: c_gammar_pp

        gammar_pp = 0.0d0
        c_gammar_pp = 0.0d0

        do i = 1, Nr5_terms
            I_val = real(Ir_r5(i), real64)
            call kahan_add(gammar_pp, c_gammar_pp, &
                           nr_r5(i) * I_val * (I_val - 1.0d0) * (pi**(Ir_r5(i) - 2)) * (tau**Jr_r5(i)))
        end do
    end function calc_gammar_pp_region5

    pure elemental function calc_gammar_t_region5(tau, pi) result(gammar_t)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        real(real64) :: gammar_t

        integer(int32) :: i
        integer(int32) :: J_val
        real(real64) :: c_gammar_t

        gammar_t = 0.0d0
        c_gammar_t = 0.0d0

        do i = 1, Nr5_terms
            J_val = real(Jr_r5(i), real64)
            call kahan_add(gammar_t, c_gammar_t, &
                           nr_r5(i) * J_val * (pi**Ir_r5(i)) * (tau**(Jr_r5(i) - 1)))
        end do
    end function calc_gammar_t_region5

    pure elemental function calc_gammar_tt_region5(tau, pi) result(gammar_tt)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Second derivative \(\gamma^r_{tt}\)
        real(real64) :: gammar_tt

        integer(int32) :: i
        integer(int32) :: J_val
        real(real64) :: c_gammar_tt

        gammar_tt = 0.0d0
        c_gammar_tt = 0.0d0

        do i = 1, Nr5_terms
            J_val = real(Jr_r5(i), real64)
            call kahan_add(gammar_tt, c_gammar_tt, &
                           nr_r5(i) * J_val * (J_val - 1.0d0) * (pi**Ir_r5(i)) * (tau**(Jr_r5(i) - 2)))
        end do
    end function calc_gammar_tt_region5

    pure elemental function calc_gammar_pt_region5(tau, pi) result(gammar_pt)
        implicit none
        !> Reduced temperature, \(\tau\)
        real(real64), intent(in) :: tau
        !> Reduced pressure, \(\pi\)
        real(real64), intent(in) :: pi
        !> Mixed derivative \(\gamma^r_{p t}\)
        real(real64) :: gammar_pt

        integer(int32) :: i
        integer(int32) :: I_val, J_val
        real(real64) :: c_gammar_pt

        gammar_pt = 0.0d0
        c_gammar_pt = 0.0d0

        do i = 1, Nr5_terms
            I_val = real(Ir_r5(i), real64)
            J_val = real(Jr_r5(i), real64)
            call kahan_add(gammar_pt, c_gammar_pt, &
                           nr_r5(i) * I_val * J_val * (pi**(Ir_r5(i) - 1)) * (tau**(Jr_r5(i) - 1)))
        end do
    end function calc_gammar_pt_region5
end submodule iapws97_region5
