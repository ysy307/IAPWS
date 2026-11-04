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

    ! Exponent range bounds for ideal-gas power table
    ! J0_r2: min=-5, max=3; derivatives need J-2 => lower bound = -7
    integer(int32), parameter :: TAU0_EXP_LO = -7, TAU0_EXP_HI = 3

    ! Exponent range bounds for residual power tables
    ! Ir_r2: min=1, max=24; derivatives need I-2 => lower bound = -1
    ! Jr_r2: min=0, max=58; derivatives need J-2 => lower bound = -2
    integer(int32), parameter :: PIR_EXP_LO = -1, PIR_EXP_HI = 24
    integer(int32), parameter :: TAUR_EXP_LO = -2, TAUR_EXP_HI = 58

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

        ! Ideal-gas and residual contributions
        real(real64) :: g0, g0_t, g0_tt
        real(real64) :: gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt

        call coef%reset()

        ! 1. Dimensionless variables (Region 2 definition)
        tau = self%T_star / T_in
        pi = P_in / self%p_star

        ! 2. Calculate ideal-gas and residual parts with pre-computed power tables
        call calc_gamma0_all_region2(tau, pi, g0, g0_t, g0_tt)
        call calc_gammar_all_region2(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt)

        ! gamma
        val_g = g0 + gr

        ! gamma_pi (ideal-gas: 1/pi)
        val_gp = 1.0d0 / pi + gr_p

        ! gamma_tau
        val_gt = g0_t + gr_t

        ! gamma_pi_pi (ideal-gas: -1/pi^2)
        val_gpp = -1.0d0 / (pi * pi) + gr_pp

        ! gamma_tau_tau
        val_gtt = g0_tt + gr_tt

        ! gamma_pi_tau (ideal-gas: 0)
        val_gpt = gr_pt

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

    !> Gibbs coefficients and cp partials share powers and coefficient traversal.
    !> Keep independent cp accumulators: their original operation order differs
    !> from the ordinary Gibbs sums and is part of numerical compatibility.
    module pure elemental subroutine calc_gamma_cp_iapws97_region2(self, T_in, P_in, coef, dcp_dT, dcp_dp)
        class(type_iapws97_region2), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        real(real64), intent(out) :: dcp_dT, dcp_dp
#if defined(__NVCOMPILER)
        ! NVHPC vectorizes the fused reductions with a different summation
        ! order. Retain bitwise compatibility with its ordinary Gibbs path.
        call self%calc_gamma(T_in,P_in,coef)
        call self%calc_cp_derivatives(T_in,P_in,dcp_dT,dcp_dp)
#else
        call evaluate_gamma_region2(self,T_in,P_in,coef,dcp_dT,dcp_dp)
#endif
    end subroutine calc_gamma_cp_iapws97_region2

    pure elemental subroutine evaluate_gamma_region2(self, T_in, P_in, coef, dcp_dT, dcp_dp)
        implicit none
        !> IAPWS-97 Region 2 object.
        class(type_iapws97_region2), intent(in) :: self
        !> temperature T, [K]
        real(real64), intent(in) :: T_in
        !> pressure P, [Pa]
        real(real64), intent(in) :: P_in
        !> IAPWS Gibbs coefficients (Dimensional)
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        real(real64), intent(out) :: dcp_dT, dcp_dp
        real(real64) :: cp_tt, cp_ttt, cp_ttp

        real(real64) :: tau, pi
        ! Local variables for dimensionless sums
        real(real64) :: val_g, val_gp, val_gt
        real(real64) :: val_gpp, val_gtt, val_gpt
        real(real64) :: RT, R_pstar

        ! Ideal-gas and residual contributions
        real(real64) :: g0, g0_t, g0_tt
        real(real64) :: gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt

        call coef%reset()

        ! 1. Dimensionless variables (Region 2 definition)
        tau = self%T_star / T_in
        pi = P_in / self%p_star

        ! 2. Calculate ideal-gas and residual parts with pre-computed power tables
        cp_ttp=0d0
        call calc_gamma0_all_cp_region2(tau, pi, g0, g0_t, g0_tt, cp_tt, cp_ttt)
        call calc_gammar_all_cp_region2(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt, cp_tt, cp_ttt, cp_ttp)

        ! gamma
        val_g = g0 + gr

        ! gamma_pi (ideal-gas: 1/pi)
        val_gp = 1.0d0 / pi + gr_p

        ! gamma_tau
        val_gt = g0_t + gr_t

        ! gamma_pi_pi (ideal-gas: -1/pi^2)
        val_gpp = -1.0d0 / (pi * pi) + gr_pp

        ! gamma_tau_tau
        val_gtt = g0_tt + gr_tt

        ! gamma_pi_tau (ideal-gas: 0)
        val_gpt = gr_pt

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

        dcp_dT=self%R*tau**2/T_in*(2d0*cp_tt+tau*cp_ttt)
        dcp_dp=-self%R*tau**2/self%p_star*cp_ttp
    end subroutine evaluate_gamma_region2

    !> Calculate all ideal-gas (gamma0) derivatives for Region 2
    !> using a pre-computed power table for tau.
    pure subroutine calc_gamma0_all_region2(tau, pi, g0, g0_t, g0_tt)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: g0, g0_t, g0_tt

        integer(int32) :: i, k
        real(real64) :: powers_tau(TAU0_EXP_LO:TAU0_EXP_HI)
        real(real64) :: inv_tau, J_val

        ! Build power table for tau: indices -7..3
        powers_tau(0) = 1.0d0
        powers_tau(1) = tau
        do k = 2, TAU0_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * tau
        end do
        inv_tau = 1.0d0 / tau
        powers_tau(-1) = inv_tau
        do k = -2, TAU0_EXP_LO, -1
            powers_tau(k) = powers_tau(k + 1) * inv_tau
        end do


        g0 = log(pi)
        g0_t = 0.0d0
        g0_tt = 0.0d0

        do i = 1, N02_terms
            J_val = real(J0_r2(i), real64)
            g0 = g0 + n0_r2(i) * powers_tau(J0_r2(i))
            g0_t = g0_t + n0_r2(i) * J_val * powers_tau(J0_r2(i) - 1)
            g0_tt = g0_tt + n0_r2(i) * J_val * (J_val - 1.0d0) * powers_tau(J0_r2(i) - 2)
        end do
    end subroutine calc_gamma0_all_region2

    pure subroutine calc_gamma0_all_cp_region2(tau, pi, g0, g0_t, g0_tt, cp_tt, cp_ttt)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: g0, g0_t, g0_tt
        real(real64), intent(inout) :: cp_tt, cp_ttt
        real(real64) :: cp_tau(-8:3), cp_term

        integer(int32) :: i, k
        real(real64) :: powers_tau(TAU0_EXP_LO:TAU0_EXP_HI)
        real(real64) :: inv_tau, J_val

        ! Build power table for tau: indices -7..3
        powers_tau(0) = 1.0d0
        powers_tau(1) = tau
        do k = 2, TAU0_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * tau
        end do
        inv_tau = 1.0d0 / tau
        powers_tau(-1) = inv_tau
        do k = -2, TAU0_EXP_LO, -1
            powers_tau(k) = powers_tau(k + 1) * inv_tau
        end do

            cp_tau(0:3)=powers_tau(0:3)
            do k=-1,-8,-1
                cp_tau(k)=cp_tau(k+1)/tau
            end do
            cp_tt=0d0; cp_ttt=0d0

        g0 = log(pi)
        g0_t = 0.0d0
        g0_tt = 0.0d0

        do i = 1, N02_terms
            J_val = real(J0_r2(i), real64)
            g0 = g0 + n0_r2(i) * powers_tau(J0_r2(i))
            g0_t = g0_t + n0_r2(i) * J_val * powers_tau(J0_r2(i) - 1)
            g0_tt = g0_tt + n0_r2(i) * J_val * (J_val - 1.0d0) * powers_tau(J0_r2(i) - 2)
                if(J0_r2(i)/=0.and.J0_r2(i)/=1) then
                    cp_term=n0_r2(i)*J_val*(J_val-1d0)
                    cp_tt=cp_tt+cp_term*cp_tau(J0_r2(i)-2)
                    if(J0_r2(i)/=2) cp_ttt=cp_ttt+cp_term*(J_val-2d0)*cp_tau(J0_r2(i)-3)
                end if
        end do
    end subroutine calc_gamma0_all_cp_region2


    !> Calculate all residual (gammar) derivatives for Region 2
    !> using pre-computed power tables for pi and (tau - 0.5).
    pure subroutine calc_gammar_all_region2(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt

        integer(int32) :: i, k
        real(real64) :: base_tau, inv_pi, inv_tau
        real(real64) :: powers_pi(PIR_EXP_LO:PIR_EXP_HI)
        real(real64) :: powers_tau(TAUR_EXP_LO:TAUR_EXP_HI)
        real(real64) :: I_val, J_val, pp, pt

        base_tau = tau - 0.5d0

        ! Build power table for pi: indices -1..24
        powers_pi(0) = 1.0d0
        powers_pi(1) = pi
        do k = 2, PIR_EXP_HI
            powers_pi(k) = powers_pi(k - 1) * pi
        end do
        inv_pi = 1.0d0 / pi
        powers_pi(-1) = inv_pi

        ! Build power table for base_tau: indices -2..58
        powers_tau(0) = 1.0d0
        powers_tau(1) = base_tau
        do k = 2, TAUR_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * base_tau
        end do
        inv_tau = 1.0d0 / base_tau
        powers_tau(-1) = inv_tau
        powers_tau(-2) = inv_tau * inv_tau

        gr = 0.0d0
        gr_p = 0.0d0
        gr_t = 0.0d0
        gr_pp = 0.0d0
        gr_tt = 0.0d0
        gr_pt = 0.0d0

        do i = 1, Nr2_terms
            I_val = real(Ir_r2(i), real64)
            J_val = real(Jr_r2(i), real64)
            pp = powers_pi(Ir_r2(i))
            pt = powers_tau(Jr_r2(i))

            gr = gr + nr_r2(i) * pp * pt

            gr_p = gr_p + &
                   nr_r2(i) * I_val * powers_pi(Ir_r2(i) - 1) * pt

            gr_t = gr_t + &
                   nr_r2(i) * pp * J_val * powers_tau(Jr_r2(i) - 1)

            gr_pp = gr_pp + &
                    nr_r2(i) * I_val * (I_val - 1.0d0) * powers_pi(Ir_r2(i) - 2) * pt

            gr_tt = gr_tt + &
                    nr_r2(i) * pp * J_val * (J_val - 1.0d0) * powers_tau(Jr_r2(i) - 2)

            gr_pt = gr_pt + &
                    nr_r2(i) * I_val * powers_pi(Ir_r2(i) - 1) * J_val * powers_tau(Jr_r2(i) - 1)
        end do
    end subroutine calc_gammar_all_region2

    pure subroutine calc_gammar_all_cp_region2(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt, cp_tt, cp_ttt, cp_ttp)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt
        real(real64), intent(inout) :: cp_tt, cp_ttt, cp_ttp
        real(real64) :: cp_term

        integer(int32) :: i, k
        real(real64) :: base_tau, inv_pi, inv_tau
        real(real64) :: powers_pi(PIR_EXP_LO:PIR_EXP_HI)
        real(real64) :: powers_tau(TAUR_EXP_LO:TAUR_EXP_HI)
        real(real64) :: I_val, J_val, pp, pt

        base_tau = tau - 0.5d0

        ! Build power table for pi: indices -1..24
        powers_pi(0) = 1.0d0
        powers_pi(1) = pi
        do k = 2, PIR_EXP_HI
            powers_pi(k) = powers_pi(k - 1) * pi
        end do
        inv_pi = 1.0d0 / pi
        powers_pi(-1) = inv_pi

        ! Build power table for base_tau: indices -2..58
        powers_tau(0) = 1.0d0
        powers_tau(1) = base_tau
        do k = 2, TAUR_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * base_tau
        end do
        inv_tau = 1.0d0 / base_tau
        powers_tau(-1) = inv_tau
        powers_tau(-2) = inv_tau * inv_tau

        gr = 0.0d0
        gr_p = 0.0d0
        gr_t = 0.0d0
        gr_pp = 0.0d0
        gr_tt = 0.0d0
        gr_pt = 0.0d0

        do i = 1, Nr2_terms
            I_val = real(Ir_r2(i), real64)
            J_val = real(Jr_r2(i), real64)
            pp = powers_pi(Ir_r2(i))
            pt = powers_tau(Jr_r2(i))

            gr = gr + nr_r2(i) * pp * pt

            gr_p = gr_p + &
                   nr_r2(i) * I_val * powers_pi(Ir_r2(i) - 1) * pt

            gr_t = gr_t + &
                   nr_r2(i) * pp * J_val * powers_tau(Jr_r2(i) - 1)

            gr_pp = gr_pp + &
                    nr_r2(i) * I_val * (I_val - 1.0d0) * powers_pi(Ir_r2(i) - 2) * pt

            gr_tt = gr_tt + &
                    nr_r2(i) * pp * J_val * (J_val - 1.0d0) * powers_tau(Jr_r2(i) - 2)

            gr_pt = gr_pt + &
                    nr_r2(i) * I_val * powers_pi(Ir_r2(i) - 1) * J_val * powers_tau(Jr_r2(i) - 1)
                if(Jr_r2(i)>=2) then
                    cp_term=nr_r2(i)*J_val*(J_val-1d0)
                    cp_tt=cp_tt+cp_term*powers_pi(Ir_r2(i))*powers_tau(Jr_r2(i)-2)
                    if(Jr_r2(i)>2) cp_ttt=cp_ttt+cp_term*(J_val-2d0)*powers_pi(Ir_r2(i))*powers_tau(Jr_r2(i)-3)
                    cp_ttp=cp_ttp+cp_term*I_val*powers_pi(Ir_r2(i)-1)*powers_tau(Jr_r2(i)-2)
                end if
        end do
    end subroutine calc_gammar_all_cp_region2


    !> Analytic cp derivatives, evaluated only when requested.
    module pure elemental subroutine calc_cp_derivatives_iapws97_region2(self, T_in, p_in, dcp_dT, dcp_dp)
        class(type_iapws97_region2), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT, dcp_dp
        real(real64) :: tau, pi, y, jval, term, gtt, gttt, gttp
        real(real64) :: xp(0:24), yp(0:58), t0p(-8:3)
        integer :: k, n, iexp, jexp

        tau = self%T_star / T_in
        pi = p_in / self%p_star
        y = tau - 0.5d0
        xp(0) = 1.d0
        do k = 1, 24
            xp(k) = xp(k-1) * pi
        end do
        yp(0) = 1.d0
        do k = 1, 58
            yp(k) = yp(k-1) * y
        end do
        t0p(0) = 1.d0
        do k = 1, 3
            t0p(k) = t0p(k-1) * tau
        end do
        do k = -1, -8, -1
            t0p(k) = t0p(k+1) / tau
        end do
        gtt = 0.d0
        gttt = 0.d0
        gttp = 0.d0
        do n = 1, N02_terms
            jexp = J0_r2(n)
            if (jexp == 0 .or. jexp == 1) cycle
            jval = real(jexp, real64)
            term = n0_r2(n) * jval * (jval-1.d0)
            gtt = gtt + term * t0p(jexp-2)
            if (jexp /= 2) gttt = gttt + term * (jval-2.d0) * t0p(jexp-3)
        end do
        do n = 1, Nr2_terms
            iexp = Ir_r2(n)
            jexp = Jr_r2(n)
            if (jexp < 2) cycle
            jval = real(jexp, real64)
            term = nr_r2(n) * jval * (jval-1.d0)
            gtt = gtt + term * xp(iexp) * yp(jexp-2)
            if (jexp > 2) gttt = gttt + term * (jval-2.d0) * xp(iexp) * yp(jexp-3)
            gttp = gttp + term * real(iexp, real64) * xp(iexp-1) * yp(jexp-2)
        end do
        ! The ideal log(pi) term has no temperature derivatives.
        dcp_dT = self%R * tau**2 / T_in * (2.d0*gtt + tau*gttt)
        dcp_dp = -self%R * tau**2 / self%p_star * gttp
    end subroutine calc_cp_derivatives_iapws97_region2

! State lanes are contiguous in each power table. The term loop remains
! ordered independently in every lane; no horizontal floating-point reduction.
    module pure recursive subroutine calc_gamma_batch_iapws97_region2(self,T_in,P_in,coef)
        use module_iapws_batch,only:IAPWS_BATCH_TILE
        class(type_iapws97_region2),intent(in)::self
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(out)::coef(:,:)
        type(type_iapws_gibbs_coefficient)::g
        integer::j
#if !defined(__NVCOMPILER)
        integer::base,width,k,term,lane
        real(real64)::xp(IAPWS_BATCH_TILE,-1:24),yp(IAPWS_BATCH_TILE,-2:58)
        real(real64)::t0(IAPWS_BATCH_TILE,-7:3),sums(IAPWS_BATCH_TILE,6),ideal(IAPWS_BATCH_TILE,3)
        real(real64)::tau(IAPWS_BATCH_TILE),pi(IAPWS_BATCH_TILE),y(IAPWS_BATCH_TILE),iy(IAPWS_BATCH_TILE)
        real(real64)::iv,jv,pp,pt,rt,r_pstar,val_g,val_gp,val_gt,val_gpp,val_gtt,val_gpt
#endif
        if(size(P_in)/=size(T_in).or.size(coef,1)/=size(T_in).or.size(coef,2)/=6) &
            error stop 'IAPWS batch: Gibbs coefficient shape'
#if defined(__NVCOMPILER)
        ! Preserve this compiler's existing within-state vector reductions.
        do j=1,size(T_in)
            call self%calc_gamma(T_in(j),P_in(j),g)
            coef(j,:)=[g%g,g%g_p,g%g_t,g%g_pp,g%g_tt,g%g_pt]
        end do
#else
        if(size(T_in)<4) then
            do j=1,size(T_in)
                call self%calc_gamma(T_in(j),P_in(j),g)
                coef(j,:)=[g%g,g%g_p,g%g_t,g%g_pp,g%g_tt,g%g_pt]
            end do
            return
        end if
        do base=1,size(T_in),IAPWS_BATCH_TILE
            width=min(IAPWS_BATCH_TILE,size(T_in)-base+1)
            !$omp simd
            do lane=1,width
                tau(lane)=self%T_star/T_in(base+lane-1);pi(lane)=P_in(base+lane-1)/self%p_star
                y(lane)=tau(lane)-0.5d0;iy(lane)=1d0/tau(lane)
                xp(lane,0)=1d0;xp(lane,1)=pi(lane);xp(lane,-1)=1d0/pi(lane)
                yp(lane,0)=1d0;yp(lane,1)=y(lane)
                yp(lane,-1)=1d0/y(lane)
                t0(lane,0)=1d0;t0(lane,1)=tau(lane);t0(lane,-1)=iy(lane)
                ideal(lane,1)=log(pi(lane));ideal(lane,2)=0d0;ideal(lane,3)=0d0
            end do
            !$omp simd
            do lane=1,width
                yp(lane,-2)=yp(lane,-1)*yp(lane,-1)
            end do
            do k=2,24
                !$omp simd
                do lane=1,width
                    xp(lane,k)=xp(lane,k-1)*pi(lane)
                end do
            end do
            do k=2,58
                !$omp simd
                do lane=1,width
                    yp(lane,k)=yp(lane,k-1)*y(lane)
                end do
            end do
            do k=2,3
                !$omp simd
                do lane=1,width
                    t0(lane,k)=t0(lane,k-1)*tau(lane)
                end do
            end do
            do k=-2,-7,-1
                !$omp simd
                do lane=1,width
                    t0(lane,k)=t0(lane,k+1)*iy(lane)
                end do
            end do
            do term=1,N02_terms
                jv=real(J0_r2(term),real64)
                !$omp simd
                do lane=1,width
                    ideal(lane,1)=ideal(lane,1)+n0_r2(term)*t0(lane,J0_r2(term))
                    ideal(lane,2)=ideal(lane,2)+n0_r2(term)*jv*t0(lane,J0_r2(term)-1)
                    ideal(lane,3)=ideal(lane,3)+n0_r2(term)*jv*(jv-1d0)*t0(lane,J0_r2(term)-2)
                end do
            end do
            sums=0d0
            do term=1,Nr2_terms
                iv=real(Ir_r2(term),real64);jv=real(Jr_r2(term),real64)
                !$omp simd private(pp,pt)
                do lane=1,width
                    pp=xp(lane,Ir_r2(term));pt=yp(lane,Jr_r2(term))
                    sums(lane,1)=sums(lane,1)+nr_r2(term)*pp*pt
                    sums(lane,2)=sums(lane,2)+nr_r2(term)*iv*xp(lane,Ir_r2(term)-1)*pt
                    sums(lane,3)=sums(lane,3)+nr_r2(term)*pp*jv*yp(lane,Jr_r2(term)-1)
                    sums(lane,4)=sums(lane,4)+nr_r2(term)*iv*(iv-1d0)*xp(lane,Ir_r2(term)-2)*pt
                    sums(lane,5)=sums(lane,5)+nr_r2(term)*pp*jv*(jv-1d0)*yp(lane,Jr_r2(term)-2)
                    sums(lane,6)=sums(lane,6)+nr_r2(term)*iv*xp(lane,Ir_r2(term)-1)*jv*yp(lane,Jr_r2(term)-1)
                end do
            end do
            r_pstar=self%R/self%p_star
            !$omp simd private(rt,val_g,val_gp,val_gt,val_gpp,val_gtt,val_gpt)
            do lane=1,width
                rt=self%R*T_in(base+lane-1)
                val_g=ideal(lane,1)+sums(lane,1)
                val_gp=1d0/pi(lane)+sums(lane,2)
                val_gt=ideal(lane,2)+sums(lane,3)
                val_gpp=-1d0/(pi(lane)*pi(lane))+sums(lane,4)
                val_gtt=ideal(lane,3)+sums(lane,5);val_gpt=sums(lane,6)
                coef(base+lane-1,1)=rt*val_g
                coef(base+lane-1,2)=(rt/self%p_star)*val_gp
                coef(base+lane-1,3)=self%R*(val_g-tau(lane)*val_gt)
                coef(base+lane-1,4)=(rt/(self%p_star**2))*val_gpp
                coef(base+lane-1,5)=(self%R*tau(lane)**2/T_in(base+lane-1))*val_gtt
                coef(base+lane-1,6)=r_pstar*(val_gp-tau(lane)*val_gpt)
            end do
        end do
#endif
    end subroutine calc_gamma_batch_iapws97_region2


end submodule iapws97_region2
