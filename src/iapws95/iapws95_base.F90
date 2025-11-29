submodule(module_iapws95) iapws95_base
    implicit none
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

    module pure elemental subroutine calc_phi_iapws95(self, tau, delta, property)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-95 helmholtz properties
        type(type_iapws_phi_property), intent(inout) :: property

        call self%calc_phi0(tau, delta, property)
        call self%calc_phir(tau, delta, property)

        property%phi = property%phi0 + property%phir
        property%phi_d = property%phi0_d + property%phir_d
        property%phi_t = property%phi0_t + property%phir_t
        property%phi_dd = property%phi0_dd + property%phir_dd
        property%phi_tt = property%phi0_tt + property%phir_tt
        property%phi_dt = property%phi0_dt + property%phir_dt
    end subroutine calc_phi_iapws95

    module pure elemental subroutine calc_phi0_iapws95(self, tau, delta, property)
        implicit none
        !> IAPWS-95 model instance
        class(type_iapws95), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced density ρ/rho_c, [-]
        real(real64), intent(in) :: delta
        !> IAPWS-95 ideal helmholtz properties
        type(type_iapws_phi_property), intent(inout) :: property

        ! Local variables
        integer(int32) :: i
        integer(int32) :: pow_i
        real(real64) :: t_val, n_val, g_val
        real(real64) :: exp_gtau, one_minus_exp, inv_one_minus_exp

        ! Kahan compensation variables
        real(real64) :: c_phi0, c_phi0_t, c_phi0_tt

        ! Initialize properties
        property%phi0 = 0.0d0
        property%phi0_t = 0.0d0
        property%phi0_tt = 0.0d0
        property%phi0_d = 0.0d0
        property%phi0_dd = 0.0d0
        property%phi0_dt = 0.0d0

        ! Initialize compensation variables
        c_phi0 = 0.0d0
        c_phi0_t = 0.0d0
        c_phi0_tt = 0.0d0

        ! Base ideal gas terms (using direct assignment for initial values)
        property%phi0 = n0_log(1) * log(delta) + n0_log(2) * log(tau)
        property%phi0_t = n0_log(2) / tau
        property%phi0_tt = -n0_log(2) / tau**2

        property%phi0_d = 1.0d0 / delta
        property%phi0_dd = -1.0d0 / delta**2
        property%phi0_dt = 0.0d0

        ! Loop 1: Power terms (Apply Kahan)
        do i = 1, size(pow)
            pow_i = pow(i)
            t_val = real(pow_i, real64)
            n_val = n0_pow(i)

            call kahan_add(property%phi0, c_phi0, n_val * (tau**t_val))

            if (pow_i /= 0) then
                call kahan_add(property%phi0_t, c_phi0_t, t_val * n_val * (tau**(pow_i - 1)))
            end if
            if (pow_i /= 0 .and. pow_i /= 1) then
                call kahan_add(property%phi0_tt, c_phi0_tt, n_val * t_val * (t_val - 1.0d0) * (tau**(pow_i - 2)))
            end if
        end do

        ! Loop 2: Exponential terms (Apply Kahan)
        do i = 1, size(n0_exp)
            n_val = n0_exp(i)
            g_val = g0(i)

            exp_gtau = exp(-g_val * tau)
            one_minus_exp = 1.0d0 - exp_gtau
            inv_one_minus_exp = 1.0d0 / one_minus_exp

            call kahan_add(property%phi0, c_phi0, n_val * log(one_minus_exp))

            call kahan_add(property%phi0_t, c_phi0_t, n_val * g_val * (inv_one_minus_exp - 1.0d0))

            call kahan_add(property%phi0_tt, c_phi0_tt, -n_val * g_val**2 * exp_gtau * (inv_one_minus_exp**2))
        end do
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
        type(type_iapws_phi_property), intent(inout) :: property

        ! --- Loop variables ---
        integer(int32) :: i

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

        ! --- Kahan Compensation Variables ---
        real(real64) :: c_phir, c_phir_d, c_phir_t, c_phir_dd, c_phir_tt, c_phir_dt

        ! --- Initialization ---
        property%phir = 0.0d0
        property%phir_d = 0.0d0
        property%phir_t = 0.0d0
        property%phir_dd = 0.0d0
        property%phir_tt = 0.0d0
        property%phir_dt = 0.0d0

        c_phir = 0.0d0
        c_phir_d = 0.0d0
        c_phir_t = 0.0d0
        c_phir_dd = 0.0d0
        c_phir_tt = 0.0d0
        c_phir_dt = 0.0d0

        ! Pre-calculate inverses and squares
        inv_delta = 1.0d0 / delta
        inv_tau = 1.0d0 / tau
        inv_delta_sq = inv_delta * inv_delta
        inv_tau_sq = inv_tau * inv_tau
        delta_sq = delta * delta
        tau_sq = tau * tau

        ! ==================================================================
        ! 1. Polynomial Terms (Terms 1-7)
        !    phi = n * delta^d * tau^t
        ! ==================================================================
        Polynomial: do i = 1, size(nr1)
            nr_val = nr1(i)
            d_val = real(d1(i), real64)
            t_val = t1(i)

            ! Base term
            bi = nr_val * (delta**d1(i)) * (tau**t_val)

            call kahan_add(property%phir, c_phir, bi)

            ! Derivatives
            call kahan_add(property%phir_d, c_phir_d, bi * (d_val * inv_delta))
            call kahan_add(property%phir_t, c_phir_t, bi * (t_val * inv_tau))

            call kahan_add(property%phir_dd, c_phir_dd, bi * (d_val * (d_val - 1.0d0) * inv_delta_sq))
            call kahan_add(property%phir_tt, c_phir_tt, bi * (t_val * (t_val - 1.0d0) * inv_tau_sq))
            call kahan_add(property%phir_dt, c_phir_dt, bi * ((d_val * t_val) * inv_delta * inv_tau))
        end do Polynomial

        ! ==================================================================
        ! 2. Exponential Terms (Terms 8-51)
        !    phi = n * delta^d * tau^t * exp(-g * delta^c)
        ! ==================================================================
        Exponential: do i = 1, size(nr2)
            nr_val = nr2(i)
            d_val = real(d2(i), real64)
            t_val = t2(i)
            c_val = real(c2(i), real64)
            g_val = gamma2(i)

            ! Common subexpressions
            delta_pow_c = delta**c2(i)
            exp_val = exp(-g_val * delta_pow_c)
            bi = nr_val * (delta**d2(i)) * (tau**t_val) * exp_val

            ! Derivative helpers
            Q = d_val - c_val * g_val * delta_pow_c
            R_val = (c_val**2) * g_val * delta_pow_c

            call kahan_add(property%phir, c_phir, bi)

            call kahan_add(property%phir_d, c_phir_d, bi * (Q * inv_delta))
            call kahan_add(property%phir_t, c_phir_t, bi * (t_val * inv_tau))

            call kahan_add(property%phir_dd, c_phir_dd, bi * ((Q**2 - Q - R_val) * inv_delta_sq))
            call kahan_add(property%phir_tt, c_phir_tt, bi * (t_val * (t_val - 1.0d0) * inv_tau_sq))
            call kahan_add(property%phir_dt, c_phir_dt, bi * (t_val * Q * inv_delta * inv_tau))
        end do Exponential

        ! ==================================================================
        ! 3. Gaussian Terms (Terms 52-54)
        !    phi = n * delta^d * tau^t * exp[ -alpha*(delta-epsilon)^2 - beta*(tau-gamma)^2 ]
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
            bi = nr_val * (delta**d3(i)) * (tau**t3(i)) * exp_val

            ! Derivative helpers
            term_d = (d_val * inv_delta) - 2.0d0 * alfa_val * delta_diff
            term_t = (t_val * inv_tau) - 2.0d0 * beta_val * tau_diff

            call kahan_add(property%phir, c_phir, bi)

            call kahan_add(property%phir_d, c_phir_d, bi * term_d)
            call kahan_add(property%phir_t, c_phir_t, bi * term_t)

            call kahan_add(property%phir_dd, c_phir_dd, bi * (term_d**2 - (d_val * inv_delta_sq) - 2.0d0 * alfa_val))
            call kahan_add(property%phir_tt, c_phir_tt, bi * (term_t**2 - (t_val * inv_tau_sq) - 2.0d0 * beta_val))
            call kahan_add(property%phir_dt, c_phir_dt, bi * (term_d * term_t))
        end do Gaussian

        ! ==================================================================
        ! 4. Non-Analytic Terms (Terms 55-56)
        !    phi = n * delta * Delta^b * Psi
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
            call kahan_add(property%phir, c_phir, bi)

            ! --- Derivatives of Components ---

            ! d(Theta)/d(delta)
            d_theta_d = AA_val * (1.0d0 / beta4_val) * diff_d * (sq_diff_d**(pow_theta - 1.0d0))

            ! d2(Theta)/d(delta)^2
            d2_theta_dd = d_theta_d / diff_d + &
                          AA_val * (1.0d0 / beta4_val) * (pow_theta - 1.0d0) * &
                          2.0d0 * sq_diff_d * (sq_diff_d**(pow_theta - 2.0d0))

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
            call kahan_add(property%phir_d, c_phir_d, &
                           (bi * inv_delta) + &
                           nr_val * delta * d_Delta_pow_b_d * psi + &
                           bi * (-2.0d0 * c_int_val * diff_d))

            ! phi_t
            call kahan_add(property%phir_t, c_phir_t, &
                           nr_val * delta * d_Delta_pow_b_t * psi + &
                           bi * (-2.0d0 * D_int_val * diff_t))

            ! phi_dd
            call kahan_add(property%phir_dd, c_phir_dd, &
                           nr_val * ( &
                           2.0d0 * d_Delta_pow_b_d * psi + &
                           delta * d2_Delta_pow_b_dd * psi + &
                           2.0d0 * delta * d_Delta_pow_b_d * d_psi_d + &
                           delta * (delta_term**b_val) * d2_psi_dd &
                           ))

            ! phi_tt
            call kahan_add(property%phir_tt, c_phir_tt, &
                           nr_val * delta * ( &
                           d2_Delta_pow_b_tt * psi + &
                           2.0d0 * d_Delta_pow_b_t * d_psi_t + &
                           (delta_term**b_val) * d2_psi_tt &
                           ))

            ! phi_dt
            call kahan_add(property%phir_dt, c_phir_dt, &
                           nr_val * ( &
                           d_Delta_pow_b_t * psi + &
                           (delta_term**b_val) * d_psi_t + &
                           delta * ( &
                           d2_Delta_pow_b_dt * psi + &
                           d_Delta_pow_b_d * d_psi_t + &
                           d_Delta_pow_b_t * d_psi_d + &
                           (delta_term**b_val) * d2_psi_dt &
                           ) &
                           ))

        end do NonAnalitic

    end subroutine calc_phir_iapws95
end submodule iapws95_base
