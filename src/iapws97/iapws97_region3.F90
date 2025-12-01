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
        type(type_iapws_phi_property), intent(inout) :: property

        call property%reset()

        ! Calculate phi and its derivatives
        property%phi = calc_phi_region3(tau, delta)
        property%phi_d = calc_phi_d_region3(tau, delta)
        property%phi_t = calc_phi_t_region3(tau, delta)
        property%phi_dd = calc_phi_dd_region3(tau, delta)
        property%phi_tt = calc_phi_tt_region3(tau, delta)
        property%phi_dt = calc_phi_dt_region3(tau, delta)
    end subroutine calc_phi_iapws97_region3

    !> Calculate phi
    !> References:
    !> - IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
    !>   Water and Steam, Table 32
    !> \(\phi = n_1 \ln(\delta) + \sum_{i=2}^{40} n_i \delta^{I_i} \tau^{J_i} \)
    pure elemental function calc_phi_region3(tau, delta) result(phi)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi

        integer(int32) :: i
        real(real64) :: c_phi

        ! Term i=1: n1 * ln(delta)
        phi = n_r3(1) * log(delta)
        c_phi = 0.0d0

        ! Terms i=2 to 40
        do i = 2, N3_terms
            call kahan_add(phi, c_phi, &
                           n_r3(i) * (delta**I_r3(i)) * (tau**J_r3(i)))
        end do
    end function calc_phi_region3

    !> Calculate phi_delta (1st derivative w.r.t delta)
    pure elemental function calc_phi_d_region3(tau, delta) result(phi_d)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi_d

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: c_phi_d

        ! i=1: n1 / delta
        phi_d = n_r3(1) / delta
        c_phi_d = 0.0d0
        do i = 2, N3_terms
            I_val = real(I_r3(i), real64)
            ! n_i * I_i * delta^(I_i-1) * tau^J_i
            call kahan_add(phi_d, c_phi_d, &
                           n_r3(i) * I_val * (delta**(I_r3(i) - 1)) * (tau**J_r3(i)))
        end do
    end function calc_phi_d_region3

    !> Calculate phi_d_delta (2nd derivative w.r.t delta)
    pure elemental function calc_phi_dd_region3(tau, delta) result(phi_dd)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi_dd

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: c_phi_dd

        ! i=1: -n1 / delta^2
        phi_dd = -n_r3(1) / (delta**2)
        c_phi_dd = 0.0d0

        do i = 2, N3_terms
            I_val = real(I_r3(i), real64)
            ! n_i * I_i * (I_i-1) * delta^(I_i-2) * tau^J_i
            call kahan_add(phi_dd, c_phi_dd, &
                           n_r3(i) * I_val * (I_val - 1.0d0) * (delta**(I_r3(i) - 2)) * (tau**J_r3(i)))
        end do
    end function calc_phi_dd_region3

    !> Calculate phi_tau (1st derivative w.r.t tau)
    pure elemental function calc_phi_t_region3(tau, delta) result(phi_t)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi_t

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_phi_t

        ! i=1 term depends only on delta, so derivative w.r.t tau is 0
        phi_t = 0.0d0
        c_phi_t = 0.0d0

        do i = 2, N3_terms
            J_val = real(J_r3(i), real64)
            ! n_i * J_i * delta^I_i * tau^(J_i-1)
            call kahan_add(phi_t, c_phi_t, &
                           n_r3(i) * J_val * (delta**I_r3(i)) * (tau**(J_r3(i) - 1)))
        end do
    end function calc_phi_t_region3

    !> Calculate phi_t_tau (2nd derivative w.r.t tau)
    pure elemental function calc_phi_tt_region3(tau, delta) result(phi_tt)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi_tt

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_phi_tt

        phi_tt = 0.0d0
        c_phi_tt = 0.0d0

        do i = 2, N3_terms
            J_val = real(J_r3(i), real64)
            ! n_i * J_i * (J_i-1) * delta^I_i * tau^(J_i-2)
            call kahan_add(phi_tt, c_phi_tt, &
                           n_r3(i) * J_val * (J_val - 1.0d0) * (delta**I_r3(i)) * (tau**(J_r3(i) - 2)))
        end do
    end function calc_phi_tt_region3

    !> Calculate phi_d_tau (Mixed derivative)
    pure elemental function calc_phi_dt_region3(tau, delta) result(phi_dt)
        implicit none
        real(real64), intent(in) :: tau
        real(real64), intent(in) :: delta
        real(real64) :: phi_dt

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: J_val
        real(real64) :: c_phi_dt

        ! i=1 term is 0 for mixed derivative
        phi_dt = 0.0d0
        c_phi_dt = 0.0d0

        do i = 2, N3_terms
            I_val = real(I_r3(i), real64)
            J_val = real(J_r3(i), real64)
            ! n_i * I_i * J_i * delta^(I_i-1) * tau^(J_i-1)
            call kahan_add(phi_dt, c_phi_dt, &
                           n_r3(i) * I_val * J_val * (delta**(I_r3(i) - 1)) * (tau**(J_r3(i) - 1)))
        end do
    end function calc_phi_dt_region3

end submodule iapws97_region3
