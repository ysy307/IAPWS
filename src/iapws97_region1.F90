submodule(module_iapws97) iapws97_region1
    !------------------------------------------------------------------------------------------
    ! Region1: Saturated liquid water (IAPWS-IF97)
    !------------------------------------------------------------------------------------------
    integer(int32), parameter :: N1_terms = 34
    integer(int32), parameter :: I_r1(N1_terms) = [ &
                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
                                 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, &
                                 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, &
                                 29, 30, 31, 32]
    integer(int32), parameter :: J_r1(N1_terms) = [ &
                                 -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, &
                                 -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, &
                                 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, &
                                 -38, -39, -40, -41]
    real(real64), parameter :: n_r1(N1_terms) = [ &
                               0.14632971213167d0, -0.84548187169114d0, -0.37563603672040d1, &
                               0.33855169168385d1, -0.95791963387872d0, 0.15772038513228d0, &
                               -0.16616417199501d-1, 0.81214629983568d-3, 0.28319080123804d-3, &
                               -0.60706301565874d-3, -0.18990068218419d-1, -0.32529748770505d-1, &
                               -0.21841717175414d-1, -0.52838357969930d-4, -0.47184321073267d-3, &
                               -0.30001780793026d-3, 0.47661393906987d-4, -0.44141845330846d-5, &
                               -0.72694996297594d-15, -0.31679644845054d-4, -0.28270797985312d-5, &
                               -0.85205128120103d-9, -0.22425281908000d-5, -0.65171222895601d-6, &
                               -0.14341729937924d-12, -0.40516996860117d-6, -0.12734301741641d-8, &
                               -0.17424871230634d-9, -0.68762131295531d-18, 0.14478307828521d-19, &
                               0.26335781662795d-22, -0.11947622640071d-22, 0.18228094581404d-23, &
                               -0.93537087292458d-25]

    ! Pre-computed exponent range bounds for power table
    ! I_r1: min=0, max=32; derivatives need I-2 => lower bound = -2
    ! J_r1: min=-41, max=17; derivatives need J-2 => lower bound = -43
    integer(int32), parameter :: PI_EXP_LO = -2, PI_EXP_HI = 32
    integer(int32), parameter :: TAU_EXP_LO = -43, TAU_EXP_HI = 17

contains
    !> Initialize the IAPWS-97 Region 1 object.
    module pure elemental subroutine initialize_type_iapws97_region1(self)
        implicit none
        !> IAPWS-97 Region 1 objects.
        class(type_iapws97_region1), intent(inout) :: self

        ! Reference values for region 1
        self%T_star = 1386.0d0 ! [K]
        self%p_star = 16.53d6 ! [Pa]
        self%R = 461.526d0 ! [J/(kg·K)]

        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_region1

    !> Calculate the Gibbs free energy \(\gamma\) for Region 1.
    module pure elemental subroutine calc_gamma_iapws97_region1(self, T_in, P_in, coef)
        implicit none
        !> IAPWS-97 Region 1 object
        class(type_iapws97_region1), intent(in) :: self
        !> Temperature T, [K]
        real(real64), intent(in) :: T_in
        !> Pressure P, [Pa]
        real(real64), intent(in) :: P_in
        !> IAPWS Gibbs coefficients
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef

        real(real64) :: tau, pi
        real(real64) :: val_g, val_gp, val_gt, val_gpp, val_gtt, val_gpt

        tau = self%T_star / T_in
        pi = P_in / self%p_star

        call coef%reset()

        ! Calculate all gamma derivatives in a single pass
        call calc_gamma_all_region1(tau, pi, val_g, val_gp, val_gt, val_gpp, val_gtt, val_gpt)

        coef%g = self%R * T_in * val_g
        coef%g_p = self%R * T_in * val_gp / self%p_star
        coef%g_t = self%R * (val_g - tau * val_gt)
        coef%g_pp = self%R * T_in * val_gpp / (self%p_star * self%p_star)
        coef%g_tt = self%R * (tau**2) * val_gtt / T_in
        coef%g_pt = self%R * (val_gp - tau * val_gpt) / self%p_star

    end subroutine calc_gamma_iapws97_region1

    !> Gibbs coefficients and cp partials share powers and coefficient traversal.
    !> Keep independent cp accumulators: their original operation order differs
    !> from the ordinary Gibbs sums and is part of numerical compatibility.
    module pure elemental subroutine calc_gamma_cp_iapws97_region1(self, T_in, P_in, coef, dcp_dT, dcp_dp)
        class(type_iapws97_region1), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        real(real64), intent(out) :: dcp_dT, dcp_dp
#if defined(__NVCOMPILER)
        ! NVHPC vectorizes the fused reductions with a different summation
        ! order. Retain bitwise compatibility with its ordinary Gibbs path.
        call self%calc_gamma(T_in,P_in,coef)
        call self%calc_cp_derivatives(T_in,P_in,dcp_dT,dcp_dp)
#else
        call evaluate_gamma_region1(self,T_in,P_in,coef,dcp_dT,dcp_dp)
#endif
    end subroutine calc_gamma_cp_iapws97_region1

    pure elemental subroutine evaluate_gamma_region1(self, T_in, P_in, coef, dcp_dT, dcp_dp)
        implicit none
        !> IAPWS-97 Region 1 object
        class(type_iapws97_region1), intent(in) :: self
        !> Temperature T, [K]
        real(real64), intent(in) :: T_in
        !> Pressure P, [Pa]
        real(real64), intent(in) :: P_in
        !> IAPWS Gibbs coefficients
        type(type_iapws_gibbs_coefficient), intent(inout) :: coef
        real(real64), intent(out) :: dcp_dT, dcp_dp
        real(real64) :: cp_tt, cp_ttt, cp_ttp

        real(real64) :: tau, pi
        real(real64) :: val_g, val_gp, val_gt, val_gpp, val_gtt, val_gpt

        tau = self%T_star / T_in
        pi = P_in / self%p_star

        call coef%reset()

        ! Calculate all gamma derivatives in a single pass
        call calc_gamma_all_cp_region1(tau, pi, val_g, val_gp, val_gt, val_gpp, val_gtt, val_gpt, cp_tt, cp_ttt, cp_ttp)

        coef%g = self%R * T_in * val_g
        coef%g_p = self%R * T_in * val_gp / self%p_star
        coef%g_t = self%R * (val_g - tau * val_gt)
        coef%g_pp = self%R * T_in * val_gpp / (self%p_star * self%p_star)
        coef%g_tt = self%R * (tau**2) * val_gtt / T_in
        coef%g_pt = self%R * (val_gp - tau * val_gpt) / self%p_star

        dcp_dT=self%R*tau**2/T_in*(2d0*cp_tt+tau*cp_ttt)
        dcp_dp=-self%R*tau**2/self%p_star*cp_ttp
    end subroutine evaluate_gamma_region1

    !> Calculate all dimensionless Gibbs free energy derivatives for Region 1
    !> in a single pass, using pre-computed power tables to eliminate
    !> repeated integer exponentiation.
    pure subroutine calc_gamma_all_region1(tau, pi, gamma, gamma_p, gamma_t, &
                                           gamma_pp, gamma_tt, gamma_pt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Resulting gamma and derivatives
        real(real64), intent(inout) :: gamma, gamma_p, gamma_t
        real(real64), intent(inout) :: gamma_pp, gamma_tt, gamma_pt

        integer(int32) :: i, k
        real(real64) :: base_pi, base_tau, inv_pi, inv_tau
        real(real64) :: powers_pi(PI_EXP_LO:PI_EXP_HI)
        real(real64) :: powers_tau(TAU_EXP_LO:TAU_EXP_HI)
        real(real64) :: I_val, J_val, pp, pt

        base_pi = 7.1d0 - pi
        base_tau = tau - 1.222d0

        ! Build power table for base_pi: indices -2..32
        powers_pi(0) = 1.0d0
        powers_pi(1) = base_pi
        do k = 2, PI_EXP_HI
            powers_pi(k) = powers_pi(k - 1) * base_pi
        end do
        inv_pi = 1.0d0 / base_pi
        powers_pi(-1) = inv_pi
        powers_pi(-2) = inv_pi * inv_pi

        ! Build power table for base_tau: indices -43..17
        powers_tau(0) = 1.0d0
        powers_tau(1) = base_tau
        do k = 2, TAU_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * base_tau
        end do
        inv_tau = 1.0d0 / base_tau
        powers_tau(-1) = inv_tau
        do k = -2, TAU_EXP_LO, -1
            powers_tau(k) = powers_tau(k + 1) * inv_tau
        end do


        ! Accumulate all 6 sums in a single loop
        gamma = 0.0d0
        gamma_p = 0.0d0
        gamma_t = 0.0d0
        gamma_pp = 0.0d0
        gamma_tt = 0.0d0
        gamma_pt = 0.0d0

        do i = 1, N1_terms
            I_val = real(I_r1(i), real64)
            J_val = real(J_r1(i), real64)
            pp = powers_pi(I_r1(i))
            pt = powers_tau(J_r1(i))

            gamma = gamma + n_r1(i) * pp * pt

            gamma_p = gamma_p + &
                      (-n_r1(i)) * I_val * powers_pi(I_r1(i) - 1) * pt

            gamma_t = gamma_t + &
                      n_r1(i) * pp * J_val * powers_tau(J_r1(i) - 1)

            gamma_pp = gamma_pp + &
                       n_r1(i) * I_val * (I_val - 1.0d0) * powers_pi(I_r1(i) - 2) * pt

            gamma_tt = gamma_tt + &
                       n_r1(i) * pp * J_val * (J_val - 1.0d0) * powers_tau(J_r1(i) - 2)

            gamma_pt = gamma_pt + &
                       (-n_r1(i)) * I_val * powers_pi(I_r1(i) - 1) * J_val * powers_tau(J_r1(i) - 1)
        end do
    end subroutine calc_gamma_all_region1

    pure subroutine calc_gamma_all_cp_region1(tau, pi, gamma, gamma_p, gamma_t, &
                                           gamma_pp, gamma_tt, gamma_pt, cp_tt, cp_ttt, cp_ttp)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Resulting gamma and derivatives
        real(real64), intent(inout) :: gamma, gamma_p, gamma_t
        real(real64), intent(inout) :: gamma_pp, gamma_tt, gamma_pt
        real(real64), intent(inout) :: cp_tt, cp_ttt, cp_ttp
        real(real64) :: cp_tau(-44:17), cp_term

        integer(int32) :: i, k
        real(real64) :: base_pi, base_tau, inv_pi, inv_tau
        real(real64) :: powers_pi(PI_EXP_LO:PI_EXP_HI)
        real(real64) :: powers_tau(TAU_EXP_LO:TAU_EXP_HI)
        real(real64) :: I_val, J_val, pp, pt

        base_pi = 7.1d0 - pi
        base_tau = tau - 1.222d0

        ! Build power table for base_pi: indices -2..32
        powers_pi(0) = 1.0d0
        powers_pi(1) = base_pi
        do k = 2, PI_EXP_HI
            powers_pi(k) = powers_pi(k - 1) * base_pi
        end do
        inv_pi = 1.0d0 / base_pi
        powers_pi(-1) = inv_pi
        powers_pi(-2) = inv_pi * inv_pi

        ! Build power table for base_tau: indices -43..17
        powers_tau(0) = 1.0d0
        powers_tau(1) = base_tau
        do k = 2, TAU_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * base_tau
        end do
        inv_tau = 1.0d0 / base_tau
        powers_tau(-1) = inv_tau
        do k = -2, TAU_EXP_LO, -1
            powers_tau(k) = powers_tau(k + 1) * inv_tau
        end do

            ! The derivative evaluator historically used division for negative
            ! powers. Keep it rather than substituting reciprocal products.
            cp_tau(0:17)=powers_tau(0:17)
            do k=-1,-44,-1
                cp_tau(k)=cp_tau(k+1)/base_tau
            end do
            cp_tt=0d0; cp_ttt=0d0; cp_ttp=0d0

        ! Accumulate all 6 sums in a single loop
        gamma = 0.0d0
        gamma_p = 0.0d0
        gamma_t = 0.0d0
        gamma_pp = 0.0d0
        gamma_tt = 0.0d0
        gamma_pt = 0.0d0

        do i = 1, N1_terms
            I_val = real(I_r1(i), real64)
            J_val = real(J_r1(i), real64)
            pp = powers_pi(I_r1(i))
            pt = powers_tau(J_r1(i))

            gamma = gamma + n_r1(i) * pp * pt

            gamma_p = gamma_p + &
                      (-n_r1(i)) * I_val * powers_pi(I_r1(i) - 1) * pt

            gamma_t = gamma_t + &
                      n_r1(i) * pp * J_val * powers_tau(J_r1(i) - 1)

            gamma_pp = gamma_pp + &
                       n_r1(i) * I_val * (I_val - 1.0d0) * powers_pi(I_r1(i) - 2) * pt

            gamma_tt = gamma_tt + &
                       n_r1(i) * pp * J_val * (J_val - 1.0d0) * powers_tau(J_r1(i) - 2)

            gamma_pt = gamma_pt + &
                       (-n_r1(i)) * I_val * powers_pi(I_r1(i) - 1) * J_val * powers_tau(J_r1(i) - 1)
                if(J_r1(i)/=0.and.J_r1(i)/=1) then
                    cp_term=n_r1(i)*J_val*(J_val-1d0)
                    cp_tt=cp_tt+cp_term*powers_pi(I_r1(i))*cp_tau(J_r1(i)-2)
                    if(J_r1(i)/=2) cp_ttt=cp_ttt+cp_term*(J_val-2d0)*powers_pi(I_r1(i))*cp_tau(J_r1(i)-3)
                    if(I_r1(i)/=0) cp_ttp=cp_ttp-cp_term*I_val*powers_pi(I_r1(i)-1)*cp_tau(J_r1(i)-2)
                end if
        end do
    end subroutine calc_gamma_all_cp_region1


    !> Analytic cp derivatives at fixed pressure and fixed temperature, respectively.
    !> This evaluator preserves the standalone cp derivative operation order.
    module pure elemental subroutine calc_cp_derivatives_iapws97_region1(self, T_in, p_in, dcp_dT, dcp_dp)
        class(type_iapws97_region1), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT, dcp_dp
        real(real64) :: tau, pi, x, y, jval, term, gtt, gttt, gttp
        real(real64) :: xp(0:32), yp(-44:17)
        integer :: k, n, iexp, jexp

        tau = self%T_star / T_in
        pi = p_in / self%p_star
        x = 7.1d0 - pi
        y = tau - 1.222d0
        xp(0) = 1.d0
        do k = 1, 32
            xp(k) = xp(k-1) * x
        end do
        yp(0) = 1.d0
        do k = 1, 17
            yp(k) = yp(k-1) * y
        end do
        do k = -1, -44, -1
            yp(k) = yp(k+1) / y
        end do
        gtt = 0.d0
        gttt = 0.d0
        gttp = 0.d0
        do n = 1, N1_terms
            iexp = I_r1(n)
            jexp = J_r1(n)
            if (jexp == 0 .or. jexp == 1) cycle
            jval = real(jexp, real64)
            term = n_r1(n) * jval * (jval-1.d0)
            gtt = gtt + term * xp(iexp) * yp(jexp-2)
            if (jexp /= 2) gttt = gttt + term * (jval-2.d0) * xp(iexp) * yp(jexp-3)
            if (iexp /= 0) gttp = gttp - term * real(iexp, real64) * xp(iexp-1) * yp(jexp-2)
        end do
        ! cp = -R*tau**2*gamma_tautau; d(tau)/dT = -tau/T.
        dcp_dT = self%R * tau**2 / T_in * (2.d0*gtt + tau*gttt)
        dcp_dp = -self%R * tau**2 / self%p_star * gttp
    end subroutine calc_cp_derivatives_iapws97_region1

! State lanes are contiguous in each power table. The term loop remains
! ordered independently in every lane; no horizontal floating-point reduction.
    module pure recursive subroutine calc_gamma_batch_iapws97_region1(self,T_in,P_in,coef)
        use module_iapws_batch,only:IAPWS_BATCH_TILE
        class(type_iapws97_region1),intent(in)::self
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(out)::coef(:,:)
        type(type_iapws_gibbs_coefficient)::g
        integer::j
#if !defined(__NVCOMPILER)
        integer::base,width,k,term,lane
        real(real64)::xp(IAPWS_BATCH_TILE,-2:32),yp(IAPWS_BATCH_TILE,-43:17)
        real(real64)::sums(IAPWS_BATCH_TILE,6),tau(IAPWS_BATCH_TILE)
        real(real64)::x(IAPWS_BATCH_TILE),y(IAPWS_BATCH_TILE),ix(IAPWS_BATCH_TILE),iy(IAPWS_BATCH_TILE)
        real(real64)::iv,jv,pp,pt
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
                tau(lane)=self%T_star/T_in(base+lane-1)
                x(lane)=7.1d0-P_in(base+lane-1)/self%p_star
                y(lane)=tau(lane)-1.222d0
                ix(lane)=1d0/x(lane);iy(lane)=1d0/y(lane)
                xp(lane,0)=1d0;xp(lane,1)=x(lane)
                xp(lane,-1)=ix(lane);xp(lane,-2)=ix(lane)*ix(lane)
                yp(lane,0)=1d0;yp(lane,1)=y(lane);yp(lane,-1)=iy(lane)
            end do
            do k=2,32
                !$omp simd
                do lane=1,width
                    xp(lane,k)=xp(lane,k-1)*x(lane)
                end do
            end do
            do k=2,17
                !$omp simd
                do lane=1,width
                    yp(lane,k)=yp(lane,k-1)*y(lane)
                end do
            end do
            do k=-2,-43,-1
                !$omp simd
                do lane=1,width
                    yp(lane,k)=yp(lane,k+1)*iy(lane)
                end do
            end do
            sums=0d0
            do term=1,N1_terms
                iv=real(I_r1(term),real64);jv=real(J_r1(term),real64)
                !$omp simd private(pp,pt)
                do lane=1,width
                    pp=xp(lane,I_r1(term));pt=yp(lane,J_r1(term))
                    sums(lane,1)=sums(lane,1)+n_r1(term)*pp*pt
                    sums(lane,2)=sums(lane,2)+(-n_r1(term))*iv*xp(lane,I_r1(term)-1)*pt
                    sums(lane,3)=sums(lane,3)+n_r1(term)*pp*jv*yp(lane,J_r1(term)-1)
                    sums(lane,4)=sums(lane,4)+n_r1(term)*iv*(iv-1d0)*xp(lane,I_r1(term)-2)*pt
                    sums(lane,5)=sums(lane,5)+n_r1(term)*pp*jv*(jv-1d0)*yp(lane,J_r1(term)-2)
                    sums(lane,6)=sums(lane,6)+(-n_r1(term))*iv*xp(lane,I_r1(term)-1)*jv*yp(lane,J_r1(term)-1)
                end do
            end do
            !$omp simd
            do lane=1,width
                coef(base+lane-1,1)=self%R*T_in(base+lane-1)*sums(lane,1)
                coef(base+lane-1,2)=self%R*T_in(base+lane-1)*sums(lane,2)/self%p_star
                coef(base+lane-1,3)=self%R*(sums(lane,1)-tau(lane)*sums(lane,3))
                coef(base+lane-1,4)=self%R*T_in(base+lane-1)*sums(lane,4)/(self%p_star*self%p_star)
                coef(base+lane-1,5)=self%R*(tau(lane)**2)*sums(lane,5)/T_in(base+lane-1)
                coef(base+lane-1,6)=self%R*(sums(lane,2)-tau(lane)*sums(lane,6))/self%p_star
            end do
        end do
#endif
    end subroutine calc_gamma_batch_iapws97_region1


end submodule iapws97_region1
