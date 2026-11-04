submodule(module_iapws97) iapws97_region3
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    integer(int32), parameter :: N3_terms = 40
    integer(int32), parameter :: I_r3(2:N3_terms) = [ &
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, &
                                 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, &
                                 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, &
                                 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
    integer(int32), parameter :: J_r3(2:N3_terms) = [ &
                                 0, 1, 2, 7, 10, 12, 23, 2, 6, &
                                 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, &
                                 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, &
                                 0, 2, 26, 2, 26, 2, 26, 0, 1, 26]
    real(real64), parameter :: n_r3(N3_terms) = [ &
                               0.10658070028513d1, -0.15732845290239d2, 0.20944396974307d2, &
                               -0.76867707878716d1, 0.26185947787954d1, -0.28080781148620d1, &
                               0.12053369696517d1, -0.84566812812502d-2, -0.12654315477714d1, &
                               -0.11524407806681d1, 0.88521043984318d0, -0.64207765181607d0, &
                               0.38493460186671d0, -0.85214708824206d0, 0.48972281541877d1, &
                               -0.30502617256965d1, 0.39420536879154d-1, 0.12558408424308d0, &
                               -0.27999329698710d0, 0.13899799569460d1, -0.20189915023570d1, &
                               -0.82147637173963d-2, -0.47596035734923d0, 0.43984074473500d-1, &
                               -0.44476435428739d0, 0.90572070719733d0, 0.70522450087967d0, &
                               0.10770512626332d0, -0.32913623258954d0, -0.50871062041158d0, &
                               -0.22175400873096d-1, 0.94260751665092d-1, 0.16436278447961d0, &
                               -0.13503372241348d-1, -0.14834345352472d-1, 0.57922953628084d-3, &
                               0.32308904703711d-2, 0.80964802996215d-4, -0.16557679795037d-3, &
                               -0.44923899061815d-4]

    ! Exponent range bounds for power tables
    ! I_r3: min=0, max=11; derivatives need I-2 => lower bound = -2
    ! J_r3: min=0, max=26; derivatives need J-2 => lower bound = -2
    integer(int32), parameter :: DELTA_EXP_LO = -2, DELTA_EXP_HI = 11
    integer(int32), parameter :: TAU_EXP_LO = -2, TAU_EXP_HI = 26

contains
    !> Pressure and density derivative for safeguarded root solvers.
    !> Only phi_d and phi_dd are needed. The specialized path preserves their
    !> recurrences, avoiding four other sums and log(delta). NVHPC instead
    !> reuses one full evaluation to retain its original vector reductions.
    module pure elemental subroutine calc_pressure_pair_iapws97_region3(self, T_in, rho_in, p, p_rho)
        class(type_iapws97_region3), intent(in) :: self
        real(real64), intent(in) :: T_in, rho_in
        real(real64), intent(out) :: p, p_rho
#if defined(__NVCOMPILER)
        type(type_iapws_helmholtz_property) :: f
        ! Its six-sum and two-sum vector reductions have different rounding.
        ! Reuse one original phi evaluation for p and p_rho rather than two;
        ! keep the root solver's pressure values and decisions bitwise stable.
        call self%calc_phi(self%T_c/T_in,rho_in/self%rho_c,f)
        call self%calc_p(T_in,rho_in,p,f)
        call self%calc_p_rho(T_in,rho_in,p_rho,f)
#else
        real(real64) :: tau, delta, inv_delta, fd, fdd, iv, term_tau
        real(real64) :: pd(DELTA_EXP_LO:DELTA_EXP_HI), pt(0:TAU_EXP_HI)
        integer :: i,k
        tau=self%T_c/T_in
        delta=rho_in/self%rho_c
        pd(0)=1d0; pd(1)=delta
        do k=2,DELTA_EXP_HI
            pd(k)=pd(k-1)*delta
        end do
        inv_delta=1d0/delta
        pd(-1)=inv_delta; pd(-2)=inv_delta*inv_delta
        pt(0)=1d0; pt(1)=tau
        do k=2,TAU_EXP_HI
            pt(k)=pt(k-1)*tau
        end do
        fd=n_r3(1)*inv_delta
        fdd=-n_r3(1)*inv_delta*inv_delta
        do i=2,N3_terms
            iv=real(I_r3(i),real64)
            term_tau=pt(J_r3(i))
            fd=fd+n_r3(i)*iv*pd(I_r3(i)-1)*term_tau
            fdd=fdd+n_r3(i)*iv*(iv-1d0)*pd(I_r3(i)-2)*term_tau
        end do
        p=rho_in*self%R*T_in*delta*fd
        p_rho=self%R*T_in*(2d0*delta*fd+delta**2*fdd)
#endif
    end subroutine calc_pressure_pair_iapws97_region3

    !> Analytic cp derivatives in (T,p), evaluated at the supplied (T,rho) state.
    !> The pressure derivative is per Pa; no density solve is performed here.
    module pure elemental subroutine calc_cp_derivatives_iapws97_region3(self, T_in, rho_in, dcp_dT, dcp_dp, prop_in)
        class(type_iapws97_region3), intent(in) :: self
        real(real64), intent(in) :: T_in, rho_in
        real(real64), intent(out) :: dcp_dT, dcp_dp
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in
        real(real64) :: tau, delta, phi, fd, ft, fdd, ftt, fdt
        real(real64) :: fddd, fddt, fdtt, fttt, a, b, a_d, a_t, b_d, b_t, c_d, c_t
        real(real64) :: pd(0:11), pt(0:26), ni, iv, jv
        integer :: i, k, id, jt

        dcp_dT = ieee_value(0.0_real64, ieee_quiet_nan)
        dcp_dp = dcp_dT
        if (.not. ieee_is_finite(T_in) .or. .not. ieee_is_finite(rho_in)) return
        if (T_in <= 0.0_real64 .or. rho_in <= 0.0_real64) return
        if (T_in == self%T_c .and. rho_in == self%rho_c) return
        tau = self%T_c / T_in
        delta = rho_in / self%rho_c
        if(present(prop_in)) then
            fd=prop_in%f_d; fdd=prop_in%f_dd
            ftt=prop_in%f_tt; fdt=prop_in%f_dt
        else
            call calc_phi_all_region3(tau, delta, phi, fd, ft, fdd, ftt, fdt)
        end if

        ! Third derivatives are evaluated only on this derivative path.
        pd(0) = 1.0_real64
        pt(0) = 1.0_real64
        do k = 1, ubound(pd, 1)
            pd(k) = pd(k - 1) * delta
        end do
        do k = 1, ubound(pt, 1)
            pt(k) = pt(k - 1) * tau
        end do
        fddd = 2.0_real64 * n_r3(1) / delta**3
        fddt = 0.0_real64
        fdtt = 0.0_real64
        fttt = 0.0_real64
        do i = 2, N3_terms
            id = I_r3(i)
            jt = J_r3(i)
            iv = real(id, real64)
            jv = real(jt, real64)
            ni = n_r3(i)
            if (id >= 3) fddd = fddd + ni * iv * (iv - 1.0_real64) * &
                (iv - 2.0_real64) * pd(id - 3) * pt(jt)
            if (id >= 2 .and. jt >= 1) fddt = fddt + ni * iv * (iv - 1.0_real64) * &
                jv * pd(id - 2) * pt(jt - 1)
            if (id >= 1 .and. jt >= 2) fdtt = fdtt + ni * iv * jv * &
                (jv - 1.0_real64) * pd(id - 1) * pt(jt - 2)
            if (jt >= 3) fttt = fttt + ni * jv * (jv - 1.0_real64) * &
                (jv - 2.0_real64) * pd(id) * pt(jt - 3)
        end do

        ! cp/R = -tau**2*phi_tt + A**2/B, p_rho = R*T*B, p_T = R*rho*A.
        a = delta * (fd - tau * fdt)
        b = 2.0_real64 * delta * fd + delta**2 * fdd
        if (.not. ieee_is_finite(b)) return
        if (b <= 0.0_real64) return
        a_t = -delta * tau * fdtt
        a_d = fd + delta * fdd - tau * fdt - delta * tau * fddt
        b_t = 2.0_real64 * delta * fdt + delta**2 * fddt
        b_d = 2.0_real64 * fd + 4.0_real64 * delta * fdd + delta**2 * fddd
        c_t = -2.0_real64 * tau * ftt - tau**2 * fttt + &
            2.0_real64 * (a / b) * a_t - (a / b)**2 * b_t
        c_d = -tau**2 * fdtt + 2.0_real64 * (a / b) * a_d - (a / b)**2 * b_d

        ! Transform cp_T|rho and cp_rho|T to cp_T|p and cp_p|T.
        dcp_dT = -(self%R / T_in) * (tau * c_t + delta * (a / b) * c_d)
        dcp_dp = c_d / (self%rho_c * T_in * b)
    end subroutine calc_cp_derivatives_iapws97_region3

    !> Initialize the IAPWS-97 Region 3 object.
    module pure elemental subroutine initialize_type_iapws97_region3(self)
        implicit none
        !> IAPWS-97 Region 3 object.
        class(type_iapws97_region3), intent(inout) :: self

        self%T_c = 647.096d0 ! [K]
        self%rho_c = 322.0d0 ! [kg/m^3]
        self%R = 461.526d0 ! [J/(kg·K)]
        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_region3

    !> Calculate the dimensionless helmholtz gibbs energy \(\phi\) for Region 3.
    module pure elemental subroutine calc_phi_iapws97_region3(self, tau, delta, property)
        implicit none
        !> IAPWS-97 Region 3 object.
        class(type_iapws97_region3), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-97 helmholtz properties
        type(type_iapws_helmholtz_property), intent(inout) :: property

        call property%reset()

        ! Calculate all phi derivatives in a single pass
        call calc_phi_all_region3(tau, delta, &
                                  property%f, property%f_d, property%f_t, &
                                  property%f_dd, property%f_tt, property%f_dt)
    end subroutine calc_phi_iapws97_region3

    !> Calculate all dimensionless Helmholtz free energy derivatives for Region 3
    !> in a single pass, using pre-computed power tables to eliminate
    !> repeated integer exponentiation.
    pure subroutine calc_phi_all_region3(tau, delta, phi, phi_d, phi_t, phi_dd, phi_tt, phi_dt)
        implicit none
        real(real64), intent(in) :: tau, delta
        real(real64), intent(inout) :: phi, phi_d, phi_t, phi_dd, phi_tt, phi_dt

        integer(int32) :: i, k
        real(real64) :: inv_delta, inv_tau
        real(real64) :: powers_delta(DELTA_EXP_LO:DELTA_EXP_HI)
        real(real64) :: powers_tau(TAU_EXP_LO:TAU_EXP_HI)
        real(real64) :: I_val, J_val, pd, pt

        ! Build power table for delta: indices -2..11
        powers_delta(0) = 1.0d0
        powers_delta(1) = delta
        do k = 2, DELTA_EXP_HI
            powers_delta(k) = powers_delta(k - 1) * delta
        end do
        inv_delta = 1.0d0 / delta
        powers_delta(-1) = inv_delta
        powers_delta(-2) = inv_delta * inv_delta

        ! Build power table for tau: indices -2..26
        powers_tau(0) = 1.0d0
        powers_tau(1) = tau
        do k = 2, TAU_EXP_HI
            powers_tau(k) = powers_tau(k - 1) * tau
        end do
        inv_tau = 1.0d0 / tau
        powers_tau(-1) = inv_tau
        powers_tau(-2) = inv_tau * inv_tau

        ! Term i=1: n1 * ln(delta)
        phi = n_r3(1) * log(delta)
        phi_d = n_r3(1) * inv_delta
        phi_t = 0.0d0
        phi_dd = -n_r3(1) * inv_delta * inv_delta
        phi_tt = 0.0d0
        phi_dt = 0.0d0

        ! Terms i=2 to 40
        do i = 2, N3_terms
            I_val = real(I_r3(i), real64)
            J_val = real(J_r3(i), real64)
            pd = powers_delta(I_r3(i))
            pt = powers_tau(J_r3(i))

            phi = phi + n_r3(i) * pd * pt

            phi_d = phi_d + &
                    n_r3(i) * I_val * powers_delta(I_r3(i) - 1) * pt

            phi_t = phi_t + &
                    n_r3(i) * J_val * pd * powers_tau(J_r3(i) - 1)

            phi_dd = phi_dd + &
                     n_r3(i) * I_val * (I_val - 1.0d0) * powers_delta(I_r3(i) - 2) * pt

            phi_tt = phi_tt + &
                     n_r3(i) * J_val * (J_val - 1.0d0) * pd * powers_tau(J_r3(i) - 2)

            phi_dt = phi_dt + &
                     n_r3(i) * I_val * J_val * powers_delta(I_r3(i) - 1) * powers_tau(J_r3(i) - 1)
        end do
    end subroutine calc_phi_all_region3

end submodule iapws97_region3
