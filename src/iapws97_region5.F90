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

    ! Exponent range bounds for ideal-gas power table
    ! J0_r5: min=-3, max=2; derivatives need J-2 => lower bound = -5
    integer(int32), parameter :: TAU0_EXP_LO = -5, TAU0_EXP_HI = 2

    ! Exponent range bounds for residual power tables
    ! Ir_r5: min=1, max=3; derivatives need I-2 => lower bound = -1
    ! Jr_r5: min=1, max=9; derivatives need J-2 => lower bound = -1
    integer(int32), parameter :: PIR_EXP_LO = -1, PIR_EXP_HI = 3
    integer(int32), parameter :: TAUR_EXP_LO = -1, TAUR_EXP_HI = 9

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

!> Calculate the dimensional Gibbs free energy g and its derivatives for Region 5.
    module pure elemental subroutine calc_gamma_iapws97_region5(self, T_in, P_in, coef)
        implicit none
        !> IAPWS-97 Region 5 object.
        class(type_iapws97_region5), intent(in) :: self
        !> temperature T, [K]
        real(real64), intent(in) :: T_in
        !> pressure P, [Pa]
        real(real64), intent(in) :: P_in
        !> IAPWS Gibbs coefficients (Dimensional output)
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

        ! 1. Dimensionless variables (Region 5 definition)
        tau = self%T_star / T_in
        pi = P_in / self%p_star

        ! 2. Calculate ideal-gas and residual parts with pre-computed power tables
        call calc_gamma0_all_region5(tau, pi, g0, g0_t, g0_tt)
        call calc_gammar_all_region5(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt)

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

        ! 3. Convert to DIMENSIONAL quantities (Adapter Logic)
        !    Same chain rules as Region 1 & 2 (due to tau = T*/T)

        RT = self%R * T_in
        R_pstar = self%R / self%p_star

        ! g [J/kg]
        coef%g = RT * val_g

        ! g_p [m^3/kg] (Specific Volume)
        coef%g_p = (RT / self%p_star) * val_gp

        ! g_t [J/(kg K)] (Negative Entropy)
        coef%g_t = self%R * (val_g - tau * val_gt)

        ! g_pp [m^3/(kg Pa)] (Derivative of volume wrt P)
        coef%g_pp = (RT / (self%p_star**2)) * val_gpp

        ! g_tt [J/(kg K^2)] (Negative Cp/T)
        coef%g_tt = (self%R * tau**2 / T_in) * val_gtt

        ! g_pt [m^3/(kg K)] (Derivative of volume wrt T)
        coef%g_pt = R_pstar * (val_gp - tau * val_gpt)

    end subroutine calc_gamma_iapws97_region5

    !> Calculate all ideal-gas (gamma0) derivatives for Region 5
    !> using a pre-computed power table for tau.
    pure subroutine calc_gamma0_all_region5(tau, pi, g0, g0_t, g0_tt)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: g0, g0_t, g0_tt

        integer(int32) :: i, k
        real(real64) :: powers_tau(TAU0_EXP_LO:TAU0_EXP_HI)
        real(real64) :: inv_tau, J_val

        ! Build power table for tau: indices -5..2
        powers_tau(0) = 1.0d0
        powers_tau(1) = tau
        powers_tau(2) = tau * tau
        inv_tau = 1.0d0 / tau
        powers_tau(-1) = inv_tau
        do k = -2, TAU0_EXP_LO, -1
            powers_tau(k) = powers_tau(k + 1) * inv_tau
        end do

        g0 = log(pi)
        g0_t = 0.0d0
        g0_tt = 0.0d0

        do i = 1, N05_terms
            J_val = real(J0_r5(i), real64)
            g0 = g0 + n0_r5(i) * powers_tau(J0_r5(i))
            g0_t = g0_t + n0_r5(i) * J_val * powers_tau(J0_r5(i) - 1)
            g0_tt = g0_tt + n0_r5(i) * J_val * (J_val - 1.0d0) * powers_tau(J0_r5(i) - 2)
        end do
    end subroutine calc_gamma0_all_region5

    !> Calculate all residual (gammar) derivatives for Region 5
    !> using pre-computed power tables for pi and tau.
    pure subroutine calc_gammar_all_region5(tau, pi, gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt)
        implicit none
        real(real64), intent(in) :: tau, pi
        real(real64), intent(inout) :: gr, gr_p, gr_t, gr_pp, gr_tt, gr_pt

        integer(int32) :: i, k
        real(real64) :: inv_pi, inv_tau
        real(real64) :: powers_pi(PIR_EXP_LO:PIR_EXP_HI)
        real(real64) :: powers_tau(TAUR_EXP_LO:TAUR_EXP_HI)
        real(real64) :: I_val, J_val, pp, pt

        ! Build power table for pi: indices -1..3
        powers_pi(0) = 1.0d0
        powers_pi(1) = pi
        do k = 2, PIR_EXP_HI
            powers_pi(k) = powers_pi(k - 1) * pi
        end do
        inv_pi = 1.0d0 / pi
        powers_pi(-1) = inv_pi

        ! Build power table for tau: indices -1..9
        powers_tau(0) = 1.0d0
        powers_tau(1) = tau
        do k = 2, TAUR_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * tau
        end do
        inv_tau = 1.0d0 / tau
        powers_tau(-1) = inv_tau

        gr = 0.0d0
        gr_p = 0.0d0
        gr_t = 0.0d0
        gr_pp = 0.0d0
        gr_tt = 0.0d0
        gr_pt = 0.0d0

        do i = 1, Nr5_terms
            I_val = real(Ir_r5(i), real64)
            J_val = real(Jr_r5(i), real64)
            pp = powers_pi(Ir_r5(i))
            pt = powers_tau(Jr_r5(i))

            gr = gr + nr_r5(i) * pp * pt

            gr_p = gr_p + &
                   nr_r5(i) * I_val * powers_pi(Ir_r5(i) - 1) * pt

            gr_t = gr_t + &
                   nr_r5(i) * J_val * pp * powers_tau(Jr_r5(i) - 1)

            gr_pp = gr_pp + &
                    nr_r5(i) * I_val * (I_val - 1.0d0) * powers_pi(Ir_r5(i) - 2) * pt

            gr_tt = gr_tt + &
                    nr_r5(i) * J_val * (J_val - 1.0d0) * pp * powers_tau(Jr_r5(i) - 2)

            gr_pt = gr_pt + &
                    nr_r5(i) * I_val * J_val * powers_pi(Ir_r5(i) - 1) * powers_tau(Jr_r5(i) - 1)
        end do
    end subroutine calc_gammar_all_region5
    !> Analytic cp derivatives, evaluated only when requested.
    module pure elemental subroutine calc_cp_derivatives_iapws97_region5(self, T_in, p_in, dcp_dT, dcp_dp)
        implicit none
        class(type_iapws97_region5), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT, dcp_dp
        real(real64) :: tau, pi, y, jval, term, gtt, gttt, gttp
        real(real64) :: xp(0:3), yp(0:9), t0p(-6:2)
        integer :: k, n, iexp, jexp

        tau = self%T_star / T_in
        pi = p_in / self%p_star
        y = tau
        xp(0) = 1.d0
        do k = 1, 3
            xp(k) = xp(k - 1) * pi
        end do
        yp(0) = 1.d0
        do k = 1, 9
            yp(k) = yp(k - 1) * y
        end do
        t0p(0) = 1.d0
        do k = 1, 2
            t0p(k) = t0p(k - 1) * tau
        end do
        do k = -1, -6, -1
            t0p(k) = t0p(k + 1) / tau
        end do
        gtt = 0.d0
        gttt = 0.d0
        gttp = 0.d0
        do n = 1, N05_terms
            jexp = J0_r5(n)
            if (jexp == 0 .or. jexp == 1) cycle
            jval = real(jexp, real64)
            term = n0_r5(n) * jval * (jval - 1.d0)
            gtt = gtt + term * t0p(jexp - 2)
            if (jexp /= 2) gttt = gttt + term * (jval - 2.d0) * t0p(jexp - 3)
        end do
        do n = 1, Nr5_terms
            iexp = Ir_r5(n)
            jexp = Jr_r5(n)
            if (jexp < 2) cycle
            jval = real(jexp, real64)
            term = nr_r5(n) * jval * (jval - 1.d0)
            gtt = gtt + term * xp(iexp) * yp(jexp - 2)
            if (jexp > 2) gttt = gttt + term * (jval - 2.d0) * xp(iexp) * yp(jexp - 3)
            gttp = gttp + term * real(iexp, real64) * xp(iexp - 1) * yp(jexp - 2)
        end do
        ! The ideal log(pi) term has no temperature derivatives.
        dcp_dT = self%R * tau**2 / T_in * (2.d0 * gtt + tau * gttt)
        dcp_dp = -self%R * tau**2 / self%p_star * gttp
    end subroutine calc_cp_derivatives_iapws97_region5

! State lanes are contiguous in each power table. The term loop remains
! ordered independently in every lane; no horizontal floating-point reduction.
    module pure recursive subroutine calc_gamma_batch_iapws97_region5(self, T_in, P_in, coef)
        use module_iapws_batch, only: IAPWS_BATCH_TILE
        implicit none
        class(type_iapws97_region5), intent(in) :: self
        real(real64), intent(in) :: T_in(:), P_in(:)
        real(real64), intent(out) :: coef(:, :)
        type(type_iapws_gibbs_coefficient) :: g
        integer :: j
#if !defined(__NVCOMPILER)
        integer :: base, width, k, term, lane
        real(real64) :: xp(IAPWS_BATCH_TILE, -1:3), yp(IAPWS_BATCH_TILE, -1:9)
        real(real64) :: t0(IAPWS_BATCH_TILE, -5:2), sums(IAPWS_BATCH_TILE, 6), ideal(IAPWS_BATCH_TILE, 3)
        real(real64) :: tau(IAPWS_BATCH_TILE), pi(IAPWS_BATCH_TILE), y(IAPWS_BATCH_TILE), iy(IAPWS_BATCH_TILE)
        real(real64) :: iv, jv, pp, pt, rt, r_pstar, val_g, val_gp, val_gt, val_gpp, val_gtt, val_gpt
#endif
        if (size(P_in) /= size(T_in) .or. size(coef, 1) /= size(T_in) .or. size(coef, 2) /= 6) &
            error stop 'IAPWS batch: Gibbs coefficient shape'
#if defined(__NVCOMPILER)
        ! Preserve this compiler's existing within-state vector reductions.
        do j = 1, size(T_in)
            call self%calc_gamma(T_in(j), P_in(j), g)
            coef(j, :) = [g%g, g%g_p, g%g_t, g%g_pp, g%g_tt, g%g_pt]
        end do
#else
        if (size(T_in) < 4) then
            do j = 1, size(T_in)
                call self%calc_gamma(T_in(j), P_in(j), g)
                coef(j, :) = [g%g, g%g_p, g%g_t, g%g_pp, g%g_tt, g%g_pt]
            end do
            return
        end if
        do base = 1, size(T_in), IAPWS_BATCH_TILE
            width = min(IAPWS_BATCH_TILE, size(T_in) - base + 1)
            !$omp simd
            do lane = 1, width
                tau(lane) = self%T_star / T_in(base + lane - 1); pi(lane) = P_in(base + lane - 1) / self%p_star
                y(lane) = tau(lane); iy(lane) = 1d0 / tau(lane)
                xp(lane, 0) = 1d0; xp(lane, 1) = pi(lane); xp(lane, -1) = 1d0 / pi(lane)
                yp(lane, 0) = 1d0; yp(lane, 1) = y(lane)
                yp(lane, -1) = 1d0 / y(lane)
                t0(lane, 0) = 1d0; t0(lane, 1) = tau(lane); t0(lane, -1) = iy(lane)
                ideal(lane, 1) = log(pi(lane)); ideal(lane, 2) = 0d0; ideal(lane, 3) = 0d0
            end do
            do k = 2, 3
                !$omp simd
                do lane = 1, width
                    xp(lane, k) = xp(lane, k - 1) * pi(lane)
                end do
            end do
            do k = 2, 9
                !$omp simd
                do lane = 1, width
                    yp(lane, k) = yp(lane, k - 1) * y(lane)
                end do
            end do
            do k = 2, 2
                !$omp simd
                do lane = 1, width
                    t0(lane, k) = t0(lane, k - 1) * tau(lane)
                end do
            end do
            do k = -2, -5, -1
                !$omp simd
                do lane = 1, width
                    t0(lane, k) = t0(lane, k + 1) * iy(lane)
                end do
            end do
            do term = 1, N05_terms
                jv = real(J0_r5(term), real64)
                !$omp simd
                do lane = 1, width
                    ideal(lane, 1) = ideal(lane, 1) + n0_r5(term) * t0(lane, J0_r5(term))
                    ideal(lane, 2) = ideal(lane, 2) + n0_r5(term) * jv * t0(lane, J0_r5(term) - 1)
                    ideal(lane, 3) = ideal(lane, 3) + n0_r5(term) * jv * (jv - 1d0) * t0(lane, J0_r5(term) - 2)
                end do
            end do
            sums = 0d0
            do term = 1, Nr5_terms
                iv = real(Ir_r5(term), real64); jv = real(Jr_r5(term), real64)
                !$omp simd private(pp,pt)
                do lane = 1, width
                    pp = xp(lane, Ir_r5(term)); pt = yp(lane, Jr_r5(term))
                    sums(lane, 1) = sums(lane, 1) + nr_r5(term) * pp * pt
                    sums(lane, 2) = sums(lane, 2) + nr_r5(term) * iv * xp(lane, Ir_r5(term) - 1) * pt
                    sums(lane, 3) = sums(lane, 3) + nr_r5(term) * jv * pp * yp(lane, Jr_r5(term) - 1)
                    sums(lane, 4) = sums(lane, 4) + nr_r5(term) * iv * (iv - 1d0) * xp(lane, Ir_r5(term) - 2) * pt
                    sums(lane, 5) = sums(lane, 5) + nr_r5(term) * jv * (jv - 1d0) * pp * yp(lane, Jr_r5(term) - 2)
                    sums(lane, 6) = sums(lane, 6) + nr_r5(term) * iv * jv * xp(lane, Ir_r5(term) - 1) * yp(lane, Jr_r5(term) - 1)
                end do
            end do
            r_pstar = self%R / self%p_star
            !$omp simd private(rt,val_g,val_gp,val_gt,val_gpp,val_gtt,val_gpt)
            do lane = 1, width
                rt = self%R * T_in(base + lane - 1)
                val_g = ideal(lane, 1) + sums(lane, 1)
                val_gp = 1d0 / pi(lane) + sums(lane, 2)
                val_gt = ideal(lane, 2) + sums(lane, 3)
                val_gpp = -1d0 / (pi(lane) * pi(lane)) + sums(lane, 4)
                val_gtt = ideal(lane, 3) + sums(lane, 5); val_gpt = sums(lane, 6)
                coef(base + lane - 1, 1) = rt * val_g
                coef(base + lane - 1, 2) = (rt / self%p_star) * val_gp
                coef(base + lane - 1, 3) = self%R * (val_g - tau(lane) * val_gt)
                coef(base + lane - 1, 4) = (rt / (self%p_star**2)) * val_gpp
                coef(base + lane - 1, 5) = (self%R * tau(lane)**2 / T_in(base + lane - 1)) * val_gtt
                coef(base + lane - 1, 6) = r_pstar * (val_gp - tau(lane) * val_gpt)
            end do
        end do
#endif
    end subroutine calc_gamma_batch_iapws97_region5

end submodule iapws97_region5
