submodule(module_iapws97) iapws97_region3
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
