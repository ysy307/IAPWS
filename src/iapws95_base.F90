submodule(module_iapws95) iapws95_base
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use module_water_phase_boundaries, only: sublimation_pressure
    implicit none

    ! Pre-computed exponent range bounds for residual power tables
    ! delta: d1 max=4, d2 max=15, d3 max=3, c2 max=6 => overall max=15
    ! tau (integer): t2 max=50, t3 max=4 => overall max=50
    integer(int32), parameter :: DELTA_POW_MAX = 15
    integer(int32), parameter :: TAU_INT_POW_MAX = 50

contains

    module pure elemental subroutine initialize_type_iapws95(self)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(inout) :: self

        self%T_c = critical_temperature
        self%rho_c = critical_density
        self%R = specific_gas_constant_water

        self%is_initialized = .true.

    end subroutine initialize_type_iapws95

    module pure elemental subroutine calc_rho_iapws95(self, T_in, P_in, rho)
        implicit none
        class(type_iapws95), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        real(real64), intent(inout) :: rho
        call density_attempt(self, T_in, P_in, rho)
    end subroutine calc_rho_iapws95

    module pure elemental subroutine solve_density_iapws95(self, T_in, P_in, rho, guess_rho)
        implicit none
        class(type_iapws95), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        real(real64), intent(out) :: rho
        real(real64), intent(in), optional :: guess_rho
        call density_attempt(self, T_in, P_in, rho, guess_rho)
        if (present(guess_rho) .and. .not. ieee_is_finite(rho)) call density_attempt(self, T_in, P_in, rho)
    end subroutine solve_density_iapws95

    pure elemental subroutine density_attempt(self, T_in, P_in, rho, guess_rho)
        implicit none
        class(type_iapws95), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        real(real64), intent(out) :: rho
        real(real64), intent(in), optional :: guess_rho
        real(real64) :: tau, theta, rt, ps, x, lo, hi, trial, f, df, best, best_error, scale, best_scale
        integer :: iteration

        rho = ieee_value(0d0, ieee_quiet_nan)
        if (.not. self%is_initialized) return
        if (.not. ieee_is_finite(T_in) .or. .not. ieee_is_finite(P_in)) return
        if (T_in < 50d0 .or. T_in > 1273d0 .or. P_in <= 0d0 .or. P_in > 1d9) return
        if (T_in < 251.165d0) then
            ps = sublimation_pressure(T_in)
            if (P_in > ps * (1d0 + 32d0 * epsilon(ps))) return
        end if
        rt = self%R * T_in
        tau = self%T_c / T_in
        lo = 0d0
        hi = 1800d0
        x = P_in / rt
        if (T_in < 251.165d0) then
            hi = 2d0 * x
        else if (tau > 1d0) then
            ! SR1-86 ancillary correlations select and seed the subcritical branch.
            ! The returned density is refined against the full Helmholtz equation.
            theta = 1d0 - 1d0 / tau
            ps = ancillary_pressure(tau, theta)
            if (T_in >= 273.16d0 .and. abs(P_in - ps) < 0.001d0 * ps) then
                call self%saturation_densities(T_in, ps, best, best_error)
                if (.not. ieee_is_finite(ps)) return
            end if
            if (P_in > ps) then
                lo = self%rho_c
                x = self%rho_c * ancillary_liquid(theta)
            else
                hi = self%rho_c
                if (P_in > 0.5d0 * ps) x = self%rho_c * ancillary_vapor(theta) * (P_in / ps)
            end if
        end if
        if (x <= lo .or. x >= hi) x = (lo + hi) * 0.5d0
        if (present(guess_rho)) then
            if (ieee_is_finite(guess_rho)) then
                ! Stay near the physical seed and avoid singular reciprocal powers.
                if (guess_rho > max(lo, 1d-120, 0.125d0 * x) .and. guess_rho < min(hi, 8d0 * x)) x = guess_rho
            end if
        end if
        best_error = huge(0d0)
        best = x
        best_scale = P_in
        do iteration = 1, 100
            call pressure_residual(x, f, df)
            if (.not. ieee_is_finite(f) .or. .not. ieee_is_finite(df)) return
            if (df > 0d0 .and. abs(f) < best_error) then
                best = x
                best_error = abs(f)
                scale = max(P_in, x * rt, x * df)
                best_scale = scale
                if (best_error <= 16d0 * epsilon(x) * scale) then
                    rho = x
                    return
                end if
            end if
            if (df <= 0d0 .and. tau > 1d0) then
                ! Move away from the unstable interval toward the selected branch.
                if (hi <= self%rho_c) then
                    hi = x
                else
                    lo = x
                end if
            else if (f > 0d0) then
                hi = x
            else
                lo = x
            end if
            trial = (lo + hi) * 0.5d0
            if (df > 0d0) then
                if (abs(f) < df * (hi - lo)) then
                    if (x - f / df > lo .and. x - f / df < hi) trial = x - f / df
                end if
            end if
            if (trial == x) exit
            x = trial
        end do
        if (best_error <= 64d0 * epsilon(best) * best_scale) rho = best

    contains
        pure subroutine pressure_residual(density, residual, slope)
            implicit none
            real(real64), intent(in) :: density
            real(real64), intent(out) :: residual, slope
            type(type_iapws_helmholtz_property) :: phi
            real(real64) :: delta, dphir
            delta = density / self%rho_c
            call self%calc_phir(tau, delta, phi)
            dphir = delta * phi%fr_d
            residual = density * rt * (1d0 + dphir) - P_in
            slope = rt * (1d0 + 2d0 * dphir + delta**2 * phi%fr_dd)
        end subroutine pressure_residual
    end subroutine density_attempt

    module pure elemental function saturation_estimate_iapws95(self, T) result(p)
        implicit none
        real(real64) :: p
        class(type_iapws95), intent(in) :: self
        real(real64), intent(in) :: T
        p = ancillary_pressure(self%T_c / T, (self%T_c - T) / self%T_c)
    end function saturation_estimate_iapws95

    pure elemental function ancillary_pressure(tau, theta) result(p)
        implicit none
        real(real64) :: p
        real(real64), intent(in) :: tau, theta
        real(real64) :: y
        y = sqrt(theta)
        p = 22.064d6 * exp(tau * theta * (-7.85951783d0 + y * 1.84408259d0 + &
            theta**2 * (-11.7866497d0 + y * 22.6807411d0 - theta * 15.9618719d0 + &
            theta**4 * y * 1.80122502d0)))
    end function ancillary_pressure

    pure elemental function ancillary_liquid(theta) result(delta)
        implicit none
        real(real64) :: delta
        real(real64), intent(in) :: theta
        real(real64) :: x
        x = theta**(1d0 / 3d0)
        delta = 1d0 + x * (1.99274064d0 + x * (1.09965342d0 - 0.510839303d0 * x**3)) &
            - 1.75493479d0 * x**16 - 45.5170352d0 * x**43 - 6.74694450d5 * x**110
    end function ancillary_liquid

    pure elemental function ancillary_vapor(theta) result(delta)
        implicit none
        real(real64) :: delta
        real(real64), intent(in) :: theta
        real(real64) :: x
        x = theta**(1d0 / 6d0)
        delta = exp(-2.03150240d0 * x**2 - 2.68302940d0 * x**4 - 5.38626492d0 * x**8 &
            - 17.2991605d0 * x**18 - 44.7586581d0 * x**37 - 63.9201063d0 * x**71)
    end function ancillary_vapor

    module pure elemental subroutine saturation_densities_iapws95(self, T, p, liquid, vapor, liquid_state, vapor_state)
        implicit none
        class(type_iapws95), intent(in) :: self
        real(real64), intent(in) :: T
        real(real64), intent(out) :: p, liquid, vapor
        type(type_iapws_property), intent(inout), optional :: liquid_state, vapor_state
        type(type_iapws_helmholtz_property) :: l, v
        real(real64) :: tau, theta, dl, dv, zl, zv, bl, bv, f1, gibbs, mismatch, ul, sl, sv, damping
        integer :: iteration

        p = ieee_value(0d0, ieee_quiet_nan)
        liquid = p
        vapor = p
        if (.not. self%is_initialized .or. .not. ieee_is_finite(T)) return
        if (T < 273.16d0 .or. T > self%T_c) return
        if (T == self%T_c) then
            p = 22.064d6
            liquid = self%rho_c
            vapor = self%rho_c
            return
        end if
        tau = self%T_c / T
        theta = (self%T_c - T) / self%T_c
        dl = ancillary_liquid(theta)
        dv = ancillary_vapor(theta)

        ! Equal reduced pressure and Gibbs energy, using the analytic Jacobian.
        do iteration = 1, 60
            call self%calc_phir(tau, dl, l)
            call self%calc_phir(tau, dv, v)
            zl = 1d0 + dl * l%fr_d
            zv = 1d0 + dv * v%fr_d
            bl = 2d0 * zl - 1d0 + dl**2 * l%fr_dd
            bv = 2d0 * zv - 1d0 + dv**2 * v%fr_dd
            f1 = dl * zl - dv * zv
            gibbs = log(dl / dv) + (l%fr - v%fr)
            mismatch = gibbs - zv * (1d0 - dv / dl)
            if (min(bl, bv) <= 0d0) return
            ! Eliminate the liquid pressure residual from the Gibbs equation.
            ! This avoids feeding its large cancellation error to the vapor step.
            sv = dv * mismatch / ((1d0 - dv / dl) * bv)
            ul = bv * sv - f1
            sl = ul / bl
            if (abs(f1) <= 128d0 * epsilon(dl) * dl * max(1d0, bl) .and. &
                abs(mismatch) <= 32d0 * epsilon(dl) * max(1d0, abs(l%fr), abs(v%fr))) then
                p = dv * gibbs / (1d0 - dv / dl) * self%rho_c * self%R * T
                liquid = dl * self%rho_c
                vapor = dv * self%rho_c
                if (present(liquid_state)) then
                    call self%calc_phi0(tau, dl, l)
                    call combine_phi(l)
                    call self%project_properties(T, liquid, liquid_state, l)
                    liquid_state%p = p
                end if
                if (present(vapor_state)) then
                    call self%calc_phi0(tau, dv, v)
                    call combine_phi(v)
                    call self%project_properties(T, vapor, vapor_state, v)
                    vapor_state%p = p
                end if
                return
            end if
            damping = 1d0
            if (sl < 0d0) damping = min(damping, 0.5d0 * (dl - 1d0) / (-sl))
            if (sv > 0d0) damping = min(damping, 0.5d0 * (1d0 - dv) / sv)
            if (sv < 0d0) damping = min(damping, 0.5d0 * dv / (-sv))
            dl = dl + damping * sl
            dv = dv + damping * sv
        end do
    end subroutine saturation_densities_iapws95

    module pure elemental subroutine calc_phi_iapws95(self, tau, delta, property)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-95 helmholtz properties
        type(type_iapws_helmholtz_property), intent(inout) :: property

        call self%calc_phi0(tau, delta, property)
        call self%calc_phir(tau, delta, property)

        call combine_phi(property)
    end subroutine calc_phi_iapws95

    pure subroutine combine_phi(property)
        implicit none
        type(type_iapws_helmholtz_property), intent(inout) :: property
        property%f = property%f0 + property%fr
        property%f_d = property%f0_d + property%fr_d
        property%f_t = property%f0_t + property%fr_t
        property%f_dd = property%f0_dd + property%fr_dd
        property%f_tt = property%f0_tt + property%fr_tt
        property%f_dt = property%f0_dt + property%fr_dt
    end subroutine combine_phi

    module pure elemental subroutine calc_phi0_iapws95(self, tau, delta, property)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-95 ideal helmholtz properties
        type(type_iapws_helmholtz_property), intent(inout) :: property

        ! Local variables
        integer(int32) :: i
        integer(int32) :: pow_i
        real(real64) :: t_val, n_val, g_val
        real(real64) :: exp_gtau, one_minus_exp, inv_one_minus_exp, q, logq
        real(real64), parameter :: ext = 0.278296458178592d0, tau_e = critical_temperature / 130d0

        ! Initialize properties
        property%f0 = 0.0d0
        property%f0_t = 0.0d0
        property%f0_tt = 0.0d0
        property%f0_d = 0.0d0
        property%f0_dd = 0.0d0
        property%f0_dt = 0.0d0

        ! Base ideal gas terms (using direct assignment for initial values)
        property%f0 = n0_log(1) * log(delta) + n0_log(2) * log(tau)
        property%f0_t = n0_log(2) / tau
        property%f0_tt = -n0_log(2) / tau**2

        property%f0_d = 1.0d0 / delta
        property%f0_dd = -1.0d0 / delta**2
        property%f0_dt = 0.0d0

        ! Loop 1: Power terms
        do i = 1, size(pow)
            pow_i = pow(i)
            t_val = real(pow_i, real64)
            n_val = n0_pow(i)

            property%f0 = property%f0 + n_val * (tau**t_val)

            if (pow_i /= 0) then
                property%f0_t = property%f0_t + t_val * n_val * (tau**(pow_i - 1))
            end if
            if (pow_i /= 0 .and. pow_i /= 1) then
                property%f0_tt = property%f0_tt + n_val * t_val * (t_val - 1.0d0) * (tau**(pow_i - 2))
            end if
        end do

        ! Loop 2: Exponential terms
        do i = 1, size(n0_exp)
            n_val = n0_exp(i)
            g_val = g0(i)

            exp_gtau = exp(-g_val * tau)
            one_minus_exp = 1.0d0 - exp_gtau
            inv_one_minus_exp = 1.0d0 / one_minus_exp

            property%f0 = property%f0 + n_val * log(one_minus_exp)

            property%f0_t = property%f0_t + n_val * g_val * (inv_one_minus_exp - 1.0d0)

            property%f0_tt = property%f0_tt + -n_val * g_val**2 * exp_gtau * (inv_one_minus_exp**2)
        end do
        ! G9-12: integrated ideal-part correction, with zero value and slope at 130 K.
        if (tau > tau_e .and. tau <= critical_temperature / 50d0) then
            q = tau / tau_e
            logq = log(q)
            property%f0 = property%f0 + ext / tau_e * &
                (0.5d0 * q**2 + 4.5d0 * (q - 1d0) - 3d0 * (q + 1d0) * logq - 0.5d0 / q)
            property%f0_t = property%f0_t + ext / tau_e**2 * &
                (q - 3d0 * logq - 3d0 / q + 0.5d0 / q**2 + 1.5d0)
            property%f0_tt = property%f0_tt + ext * (1d0 / tau_e - 1d0 / tau)**3
        end if
    end subroutine calc_phi0_iapws95

    module pure elemental subroutine calc_phir_iapws95(self, tau, delta, property)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-95 residual helmholtz properties (phi and derivatives)
        type(type_iapws_helmholtz_property), intent(inout) :: property

        ! --- Loop variables ---
        integer(int32) :: i, k

        ! --- Pre-computed power tables ---
        real(real64) :: powers_delta(0:DELTA_POW_MAX)
        real(real64) :: powers_tau_int(0:TAU_INT_POW_MAX)

        ! --- Helper variables for optimization ---
        real(real64) :: inv_delta, inv_tau, inv_delta_sq, inv_tau_sq
        real(real64) :: delta_sq, tau_sq
        real(real64) :: bi, term_d, term_t

        ! --- Polynomial & Exponential terms vars ---
        real(real64) :: d_val, t_val, nr_val, c_val, g_val
        real(real64) :: delta_pow_c, exp_val
        real(real64) :: Q, R_val

        ! --- Gaussian terms vars ---
        real(real64) :: alfa_val, beta_val, gamma_val, eps_val
        real(real64) :: delta_diff, tau_diff, sq_delta_diff, sq_tau_diff
        real(real64) :: exp_arg

        ! --- Non-Analytic terms vars ---
        real(real64) :: a_val, b_val, BB_val, AA_val, beta4_val
        real(real64) :: c_int_val, D_int_val
        real(real64) :: diff_d, diff_t, sq_diff_d, sq_diff_t
        real(real64) :: theta, delta_term, psi
        real(real64) :: pow_theta
        real(real64) :: d_theta_d, d_psi_d, d_psi_t
        real(real64) :: d2_theta_dd, d2_psi_dd, d2_psi_tt, d2_psi_dt
        real(real64) :: d_Delta_d, d_Delta_t, d2_Delta_dd, d2_Delta_tt, d2_Delta_dt
        real(real64) :: d_Delta_pow_b_d, d_Delta_pow_b_t
        real(real64) :: d2_Delta_pow_b_dd, d2_Delta_pow_b_tt, d2_Delta_pow_b_dt

        ! --- Initialization ---
        property%fr = 0.0d0
        property%fr_d = 0.0d0
        property%fr_t = 0.0d0
        property%fr_dd = 0.0d0
        property%fr_tt = 0.0d0
        property%fr_dt = 0.0d0

        ! Pre-calculate inverses and squares
        inv_delta = 1.0d0 / delta
        inv_tau = 1.0d0 / tau
        inv_delta_sq = inv_delta * inv_delta
        inv_tau_sq = inv_tau * inv_tau
        delta_sq = delta * delta
        tau_sq = tau * tau

        ! Build power table for delta: indices 0..15
        powers_delta(0) = 1.0d0
        powers_delta(1) = delta
        do k = 2, DELTA_POW_MAX
            powers_delta(k) = powers_delta(k - 1) * delta
        end do

        ! Build power table for tau (integer exponents): indices 0..50
        powers_tau_int(0) = 1.0d0
        powers_tau_int(1) = tau
        do k = 2, TAU_INT_POW_MAX
            powers_tau_int(k) = powers_tau_int(k - 1) * tau
        end do

        ! ==================================================================
        ! 1. Polynomial Terms (Terms 1-7)
        !    phi = n * delta^d * tau^t
        !    Note: t1 contains non-integer exponents, so tau**t_val
        !    cannot use the integer power table.
        ! ==================================================================
        Polynomial: do i = 1, size(nr1)
            nr_val = nr1(i)
            d_val = real(d1(i), real64)
            t_val = t1(i)

            ! Base term: use power table for delta, keep tau**t_val as-is
            bi = nr_val * powers_delta(d1(i)) * (tau**t_val)

            property%fr = property%fr + bi

            ! Derivatives
            property%fr_d = property%fr_d + bi * (d_val * inv_delta)
            property%fr_t = property%fr_t + bi * (t_val * inv_tau)

            property%fr_dd = property%fr_dd + bi * (d_val * (d_val - 1.0d0) * inv_delta_sq)
            property%fr_tt = property%fr_tt + bi * (t_val * (t_val - 1.0d0) * inv_tau_sq)
            property%fr_dt = property%fr_dt + bi * (d_val * t_val * inv_delta * inv_tau)
        end do Polynomial

        ! ==================================================================
        ! 2. Exponential Terms (Terms 8-51)
        !    phi = n * delta^d * tau^t * exp(-g * delta^c)
        !    All exponents d2, t2, c2 are integer => use power tables.
        ! ==================================================================
        Exponential: do i = 1, size(nr2)
            nr_val = nr2(i)
            d_val = real(d2(i), real64)
            t_val = real(t2(i), real64)
            c_val = real(c2(i), real64)
            g_val = real(gamma2(i), real64)

            ! Common subexpressions using power tables
            delta_pow_c = powers_delta(c2(i))
            exp_val = exp(-g_val * delta_pow_c)
            bi = nr_val * powers_delta(d2(i)) * powers_tau_int(t2(i)) * exp_val

            ! Derivative helpers
            Q = d_val - c_val * g_val * delta_pow_c
            R_val = (c_val**2) * g_val * delta_pow_c

            property%fr = property%fr + bi

            property%fr_d = property%fr_d + bi * (Q * inv_delta)
            property%fr_t = property%fr_t + bi * (t_val * inv_tau)

            property%fr_dd = property%fr_dd + bi * ((Q**2 - Q - R_val) * inv_delta_sq)
            property%fr_tt = property%fr_tt + bi * (t_val * (t_val - 1.0d0) * inv_tau_sq)
            property%fr_dt = property%fr_dt + bi * (t_val * Q * inv_delta * inv_tau)
        end do Exponential

        ! ==================================================================
        ! 3. Gaussian Terms (Terms 52-54)
        !    phi = n * delta^d * tau^t * exp[ -alpha*(delta-epsilon)^2 - beta*(tau-gamma)^2 ]
        !    d3 and t3 are integer => use power tables.
        ! ==================================================================
        Gaussian: do i = 1, size(nr3)
            nr_val = nr3(i)
            d_val = real(d3(i), real64)
            t_val = real(t3(i), real64)
            alfa_val = alfa3(i)
            beta_val = beta3(i)
            gamma_val = gamma3(i)
            eps_val = epsilon3(i)

            delta_diff = delta - eps_val
            tau_diff = tau - gamma_val
            sq_delta_diff = delta_diff * delta_diff
            sq_tau_diff = tau_diff * tau_diff

            exp_arg = -alfa_val * sq_delta_diff - beta_val * sq_tau_diff
            exp_val = exp(exp_arg)
            bi = nr_val * powers_delta(d3(i)) * powers_tau_int(t3(i)) * exp_val

            ! Derivative helpers
            term_d = (d_val * inv_delta) - 2.0d0 * alfa_val * delta_diff
            term_t = (t_val * inv_tau) - 2.0d0 * beta_val * tau_diff

            property%fr = property%fr + bi

            property%fr_d = property%fr_d + bi * term_d
            property%fr_t = property%fr_t + bi * term_t

            property%fr_dd = property%fr_dd + bi * (term_d**2 - (d_val * inv_delta_sq) - 2.0d0 * alfa_val)
            property%fr_tt = property%fr_tt + bi * (term_t**2 - (t_val * inv_tau_sq) - 2.0d0 * beta_val)
            property%fr_dt = property%fr_dt + bi * (term_d * term_t)
        end do Gaussian

        ! ==================================================================
        ! 4. Non-Analytic Terms (Terms 55-56)
        !    phi = n * delta * Delta^b * Psi
        !    Uses real-valued exponents => cannot use integer power tables.
        ! ==================================================================
        NonAnalitic: do i = 1, size(nr4)
            nr_val = nr4(i)
            a_val = a4(i)
            b_val = b4(i)
            BB_val = B(i)
            c_int_val = real(C(i), real64)
            D_int_val = real(D(i), real64)
            AA_val = A(i)
            beta4_val = beta4(i)

            diff_d = delta - 1.0d0
            diff_t = tau - 1.0d0
            sq_diff_d = diff_d * diff_d
            sq_diff_t = diff_t * diff_t

            ! Theta
            pow_theta = 1.0d0 / (2.0d0 * beta4_val)
            theta = (1.0d0 - tau) + AA_val * (sq_diff_d**pow_theta)

            ! Delta (delta_term)
            delta_term = theta**2 + BB_val * (sq_diff_d**a_val)
            ! Psi
            psi = exp(-c_int_val * sq_diff_d - D_int_val * sq_diff_t)

            ! Base phi
            bi = nr_val * delta * (delta_term**b_val) * psi
            property%fr = property%fr + bi

            ! --- Derivatives of Components ---

            ! d(Theta)/d(delta)
            d_theta_d = AA_val * (1.0d0 / beta4_val) * diff_d * (sq_diff_d**(pow_theta - 1.0d0))

            ! d2(Theta)/d(delta)^2
            ! This combined power has a finite zero limit at the critical density.
            d2_theta_dd = AA_val / beta4_val * (2.0d0 * pow_theta - 1.0d0) * &
                          sq_diff_d**(pow_theta - 1.0d0)

            ! d(Psi)/d(delta), d(Psi)/d(tau)
            d_psi_d = psi * (-2.0d0 * c_int_val * diff_d)
            d_psi_t = psi * (-2.0d0 * D_int_val * diff_t)

            ! d2(Psi)
            d2_psi_dd = d_psi_d * (-2.0d0 * c_int_val * diff_d) + psi * (-2.0d0 * c_int_val)
            d2_psi_tt = d_psi_t * (-2.0d0 * D_int_val * diff_t) + psi * (-2.0d0 * D_int_val)
            d2_psi_dt = d_psi_d * (-2.0d0 * D_int_val * diff_t)

            ! d(Delta)
            d_Delta_d = 2.0d0 * theta * d_theta_d + &
                        2.0d0 * BB_val * a_val * diff_d * (sq_diff_d**(a_val - 1.0d0))
            ! d(Delta)/dt = 2*Theta*(-1)
            d_Delta_t = -2.0d0 * theta

            ! d2(Delta)
            d2_Delta_dd = 2.0d0 * (d_theta_d**2 + theta * d2_theta_dd) + &
                          2.0d0 * BB_val * a_val * (sq_diff_d**(a_val - 1.0d0)) + &
                          4.0d0 * BB_val * a_val * (a_val - 1.0d0) * sq_diff_d * (sq_diff_d**(a_val - 2.0d0))
            d2_Delta_tt = 2.0d0 ! (-1 * -1 * 2)
            d2_Delta_dt = -2.0d0 * d_theta_d

            ! d(Delta^b)
            if (delta_term > 0.0d0) then
                d_Delta_pow_b_d = b_val * (delta_term**(b_val - 1.0d0)) * d_Delta_d
                d_Delta_pow_b_t = b_val * (delta_term**(b_val - 1.0d0)) * d_Delta_t

                ! d2(Delta^b) terms
                d2_Delta_pow_b_dd = b_val * ((b_val - 1.0d0) * (delta_term**(b_val - 2.0d0)) * d_Delta_d**2 + &
                                             (delta_term**(b_val - 1.0d0)) * d2_Delta_dd)
                d2_Delta_pow_b_tt = b_val * ((b_val - 1.0d0) * (delta_term**(b_val - 2.0d0)) * d_Delta_t**2 + &
                                             (delta_term**(b_val - 1.0d0)) * d2_Delta_tt)
                d2_Delta_pow_b_dt = b_val * ((b_val - 1.0d0) * (delta_term**(b_val - 2.0d0)) * d_Delta_d * d_Delta_t + &
                                             (delta_term**(b_val - 1.0d0)) * d2_Delta_dt)
            else
                d_Delta_pow_b_d = 0.0d0
                d_Delta_pow_b_t = 0.0d0
                d2_Delta_pow_b_dd = 0.0d0
                d2_Delta_pow_b_tt = 0.0d0
                d2_Delta_pow_b_dt = 0.0d0
            end if

            ! --- Assembly of Derivatives ---

            ! phi_d
            property%fr_d = property%fr_d + &
                            (bi * inv_delta) + &
                            nr_val * delta * d_Delta_pow_b_d * psi + &
                            bi * (-2.0d0 * c_int_val * diff_d)

            ! phi_t
            property%fr_t = property%fr_t + &
                            nr_val * delta * d_Delta_pow_b_t * psi + &
                            bi * (-2.0d0 * D_int_val * diff_t)

            ! phi_dd
            property%fr_dd = property%fr_dd + &
                             nr_val * ( &
                             2.0d0 * d_Delta_pow_b_d * psi + &
                             2.0d0 * (delta_term**b_val) * d_psi_d + &
                             delta * d2_Delta_pow_b_dd * psi + &
                             2.0d0 * delta * d_Delta_pow_b_d * d_psi_d + &
                             delta * (delta_term**b_val) * d2_psi_dd &
                             )

            ! phi_tt
            property%fr_tt = property%fr_tt + &
                             nr_val * delta * ( &
                             d2_Delta_pow_b_tt * psi + &
                             2.0d0 * d_Delta_pow_b_t * d_psi_t + &
                             (delta_term**b_val) * d2_psi_tt &
                             )

            ! phi_dt
            property%fr_dt = property%fr_dt + &
                             nr_val * ( &
                             d_Delta_pow_b_t * psi + &
                             (delta_term**b_val) * d_psi_t + &
                             delta * ( &
                             d2_Delta_pow_b_dt * psi + &
                             d_Delta_pow_b_d * d_psi_t + &
                             d_Delta_pow_b_t * d_psi_d + &
                             (delta_term**b_val) * d2_psi_dt &
                             ) &
                             )

        end do NonAnalitic

    end subroutine calc_phir_iapws95

end submodule iapws95_base
