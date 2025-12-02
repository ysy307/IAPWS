submodule(module_iapws97) iapws97_region2
    implicit none
    !------------------------------------------------------------------------------------------
    ! Region2: Superheated Steam (IAPWS-IF97)
    !------------------------------------------------------------------------------------------
    integer(int32), parameter :: N02_terms = 9
    integer(int32), parameter :: J0_r2(N02_terms) = [0, 1, -5, -4, -3, -2, -1, 2, 3]
    real(real64), parameter :: n0_r2(N02_terms) = &
                               [-0.96927686500217d+1, 0.10086655968018d+2, -0.56087911283020d-2, &
                                0.71452738081455d-1, -0.40710498223928d+0, 0.14240819171444d+1, &
                                -0.43839511319450d+1, -0.28408632460772d+0, 0.21268463753307d-1]
    integer(int32), parameter :: Nr2_terms = 43
    integer(int32), parameter :: Ir_r2(Nr2_terms) = [ &
                                 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, &
                                 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, &
                                 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, &
                                 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, &
                                 24, 24, 24]
    integer(int32), parameter :: Jr_r2(Nr2_terms) = [ &
                                 0, 1, 2, 3, 6, 1, 2, 4, 7, 36, &
                                 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, &
                                 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, &
                                 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, &
                                 26, 40, 58]
    real(real64), parameter :: nr_r2(Nr2_terms) = [ &
                               -0.17731742473213d-02, -0.17834862292358d-01, -0.45996013696365d-01, &
                               -0.57581259083432d-01, -0.50325278727930d-01, -0.33032641670203d-04, &
                               -0.18948987516315d-03, -0.39392777243355d-02, -0.43797295650573d-01, &
                               -0.26674547914087d-04, 0.20481737692309d-07, 0.43870667284435d-06, &
                               -0.32277677238570d-04, -0.15033924542148d-02, -0.40668253562649d-01, &
                               -0.78847309559367d-09, 0.12790717852285d-07, 0.48225372718507d-06, &
                               0.22922076337661d-05, -0.16714766451061d-10, -0.21171472321355d-02, &
                               -0.23895741934104d+02, -0.59059564324270d-17, -0.12621808899101d-05, &
                               -0.38946842435739d-01, 0.11256211360459d-10, -0.82311340897998d+01, &
                               0.19809712802088d-07, 0.10406965210174d-18, -0.10234747095929d-12, &
                               -0.10018179379511d-08, -0.80882908646985d-10, 0.10693031879409d+00, &
                               -0.33662250574171d+00, 0.89185845355421d-24, 0.30629316876232d-12, &
                               -0.42002467698208d-05, -0.59056029685639d-25, 0.37826947613457d-05, &
                               -0.12768608934681d-14, 0.73087610595061d-28, 0.55414715350778d-16, &
                               -0.94369707241210d-06]
contains

    !> Initialize the IAPWS-97 Region 2 object.
    module pure elemental subroutine initialize_type_iapws97_region2(self)
        implicit none
        !> IAPWS-97 Region 2 object.
        class(type_iapws97_region2), intent(inout) :: self

        self%T_star = 540.0d0 ! [K]
        self%p_star = 1.0d6 ! [Pa]
        self%R = 461.526d0 ! [J/(kg·K)]

        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_region2

!> Calculate the dimensional Gibbs free energy g and its derivatives for Region 2.
    module pure elemental subroutine calc_gamma_iapws97_region2(self, T_in, P_in, coef)
        implicit none
        !> IAPWS-97 Region 2 object.
        class(type_iapws97_region2), intent(in) :: self
        !> temperature T, [K]
        real(real64), intent(in) :: T_in
        !> pressure P, [Pa]
        real(real64), intent(in) :: P_in
        !> IAPWS Gibbs coefficients (Dimensional)
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef

        real(real64) :: tau, pi
        ! Local variables for dimensionless sums
        real(real64) :: val_g, val_gp, val_gt
        real(real64) :: val_gpp, val_gtt, val_gpt
        real(real64) :: RT, R_pstar

        call coef%reset()

        ! 1. Dimensionless variables (Region 2 definition)
        tau = self%T_star / T_in
        pi = P_in / self%p_star

        ! 2. Summing Ideal-gas part (0) and Residual part (r)
        !    Region 2: gamma = gamma0 + gammar

        ! gamma
        val_g = calc_gamma0_region2(tau, pi) + calc_gammar_region2(tau, pi)

        ! gamma_pi
        val_gp = calc_gamma0_p_region2(tau, pi) + calc_gammar_p_region2(tau, pi)

        ! gamma_tau
        val_gt = calc_gamma0_t_region2(tau, pi) + calc_gammar_t_region2(tau, pi)

        ! gamma_pi_pi
        val_gpp = calc_gamma0_pp_region2(tau, pi) + calc_gammar_pp_region2(tau, pi)

        ! gamma_tau_tau
        val_gtt = calc_gamma0_tt_region2(tau, pi) + calc_gammar_tt_region2(tau, pi)

        ! gamma_pi_tau
        val_gpt = calc_gamma0_pt_region2(tau, pi) + calc_gammar_pt_region2(tau, pi)

        ! 3. Convert to DIMENSIONAL quantities (The Adapter Logic)
        !    Using chain rules for tau = T*/T and pi = P/p*

        RT = self%R * T_in
        R_pstar = self%R / self%p_star ! Pre-calc R/p*

        ! g [J/kg] = R * T * gamma
        coef%g = RT * val_g

        ! g_p [m^3/kg] = (RT / p*) * gamma_pi
        coef%g_p = (RT / self%p_star) * val_gp

        ! g_t [J/(kg K)] = R * (gamma - tau * gamma_tau)  <-- Entropy term (-s)
        coef%g_t = self%R * (val_g - tau * val_gt)

        ! g_pp [m^3/(kg Pa)] = (RT / p*^2) * gamma_pi_pi
        coef%g_pp = (RT / (self%p_star**2)) * val_gpp

        ! g_tt [J/(kg K^2)] = (R * tau^2 / T) * gamma_tau_tau <-- Relates to -Cp/T
        coef%g_tt = (self%R * tau**2 / T_in) * val_gtt

        ! g_pt [m^3/(kg K)] = (R / p*) * (gamma_pi - tau * gamma_pi_tau)
        coef%g_pt = R_pstar * (val_gp - tau * val_gpt)

    end subroutine calc_gamma_iapws97_region2

    pure elemental function calc_gamma0_region2(tau, pi) result(gamma0)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Resulting \(\gamma^o\) value
        real(real64) :: gamma0

        integer(int32) :: i
        real(real64) :: c_gamma0

        gamma0 = log(pi)
        c_gamma0 = 0.0d0

        do i = 1, N02_terms
            call kahan_add(gamma0, c_gamma0, &
                           n0_r2(i) * tau**J0_r2(i))
        end do

    end function calc_gamma0_region2

    pure elemental function calc_gamma0_p_region2(tau, pi) result(gamma0_p)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^o_{\pi}\)
        real(real64) :: gamma0_p

        gamma0_p = 1.0d0 / pi

    end function calc_gamma0_p_region2

    pure elemental function calc_gamma0_pp_region2(tau, pi) result(gamma0_pp)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^o_{\pi}\)
        real(real64) :: gamma0_pp

        gamma0_pp = -1.0d0 / (pi * pi)

    end function calc_gamma0_pp_region2

    pure elemental function calc_gamma0_t_region2(tau, pi) result(gamma0_t)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi}\)
        real(real64) :: gamma0_t

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma0_t

        gamma0_t = 0.0d0
        c_gamma0_t = 0.0d0

        do i = 1, N02_terms
            J_val = real(J0_r2(i), real64)
            call kahan_add(gamma0_t, c_gamma0_t, &
                           n0_r2(i) * J_val * tau**(J0_r2(i) - 1))
        end do

    end function calc_gamma0_t_region2

    pure elemental function calc_gamma0_tt_region2(tau, pi) result(gamma0_tt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi}\)
        real(real64) :: gamma0_tt

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma_tt

        gamma0_tt = 0.0d0
        c_gamma_tt = 0.0d0

        do i = 1, N02_terms
            J_val = real(J0_r2(i), real64)
            call kahan_add(gamma0_tt, c_gamma_tt, &
                           n0_r2(i) * J_val * (J_val - 1.0d0) * tau**(J0_r2(i) - 2))
        end do

    end function calc_gamma0_tt_region2

    pure elemental function calc_gamma0_pt_region2(tau, pi) result(gamma0_pt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi}\)
        real(real64) :: gamma0_pt

        gamma0_pt = 0.0d0

    end function calc_gamma0_pt_region2

    pure elemental function calc_gammar_region2(tau, pi) result(gammar)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Resulting \(\gamma^r\) value
        real(real64) :: gammar

        integer(int32) :: i
        real(real64) :: c_gammar

        gammar = 0.0d0
        c_gammar = 0.0d0

        do i = 1, Nr2_terms
            call kahan_add(gammar, c_gammar, &
                           nr_r2(i) * pi**Ir_r2(i) * (tau - 0.5d0)**Jr_r2(i))
        end do

    end function calc_gammar_region2

    pure elemental function calc_gammar_p_region2(tau, pi) result(gammar_p)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi}\)
        real(real64) :: gammar_p

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: c_gammar_p

        gammar_p = 0.0d0
        c_gammar_p = 0.0d0

        do i = 1, Nr2_terms
            I_val = real(Ir_r2(i), real64)
            call kahan_add(gammar_p, c_gammar_p, &
                           nr_r2(i) * I_val * pi**(Ir_r2(i) - 1) * (tau - 0.5d0)**Jr_r2(i))
        end do

    end function calc_gammar_p_region2

    pure elemental function calc_gammar_pp_region2(tau, pi) result(gammar_pp)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi\pi}\)
        real(real64) :: gammar_pp

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: c_gammar_pp

        gammar_pp = 0.0d0
        c_gammar_pp = 0.0d0

        do i = 1, Nr2_terms
            I_val = real(Ir_r2(i), real64)
            call kahan_add(gammar_pp, c_gammar_pp, &
                           nr_r2(i) * I_val * (I_val - 1.0d0) * pi**(Ir_r2(i) - 2) * (tau - 0.5d0)**Jr_r2(i))
        end do

    end function calc_gammar_pp_region2

    pure elemental function calc_gammar_t_region2(tau, pi) result(gammar_t)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\tau}\)
        real(real64) :: gammar_t

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gammar_t

        gammar_t = 0.0d0
        c_gammar_t = 0.0d0

        do i = 1, Nr2_terms
            J_val = real(Jr_r2(i), real64)
            call kahan_add(gammar_t, c_gammar_t, &
                           nr_r2(i) * pi**Ir_r2(i) * J_val * (tau - 0.5d0)**(Jr_r2(i) - 1))
        end do

    end function calc_gammar_t_region2

    pure elemental function calc_gammar_tt_region2(tau, pi) result(gammar_tt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\tau\tau}\)
        real(real64) :: gammar_tt

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gammar_tt

        gammar_tt = 0.0d0
        c_gammar_tt = 0.0d0

        do i = 1, Nr2_terms
            J_val = real(Jr_r2(i), real64)
            call kahan_add(gammar_tt, c_gammar_tt, &
                           nr_r2(i) * pi**Ir_r2(i) * J_val * (J_val - 1.0d0) * (tau - 0.5d0)**(Jr_r2(i) - 2))
        end do

    end function calc_gammar_tt_region2

    pure elemental function calc_gammar_pt_region2(tau, pi) result(gammar_pt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma^r_{\pi\tau}\)
        real(real64) :: gammar_pt

        integer(int32) :: i
        real(real64) :: I_val, J_val
        real(real64) :: c_gammar_pt

        gammar_pt = 0.0d0
        c_gammar_pt = 0.0d0

        do i = 1, Nr2_terms
            I_val = real(Ir_r2(i), real64)
            J_val = real(Jr_r2(i), real64)
            call kahan_add(gammar_pt, c_gammar_pt, &
                           nr_r2(i) * I_val * pi**(Ir_r2(i) - 1) * J_val * (tau - 0.5d0)**(Jr_r2(i) - 1))
        end do

    end function calc_gammar_pt_region2

end submodule iapws97_region2
