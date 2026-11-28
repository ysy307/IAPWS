module iapws
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use :: module_iapws, only: type_iapws_property, IAPWS_BATCH_OK, IAPWS_BATCH_INVALID
    use :: module_iapws95, only: type_iapws95
    use :: module_iapws97, only: type_iapws97
    use :: module_iapws06_iceIh, only: type_iapws06
    use module_supercooled_water, only: in_supercooled_domain, evaluate_supercooled
    use module_water_transport, only: viscosity, thermal_conductivity, surface_tension
    use module_water_phase_boundaries, only: melting_pressure, sublimation_pressure, triple_T, triple_p, &
                                            ice_ih_melting_slope, sublimation_slope, ICE_IH, ICE_III, ICE_V, ICE_VI, ICE_VII
    implicit none
    private

    public :: IAPWS_BATCH_OK, IAPWS_BATCH_INVALID, type_iapws_property
    public :: type_iapws95, type_iapws97, type_iapws06
    public :: new_water, melting_pressure, sublimation_pressure, ICE_IH, ICE_III, ICE_V, ICE_VI, ICE_VII
    public :: viscosity, thermal_conductivity, surface_tension

    integer, parameter, public :: IF97 = 97, IAPWS95 = 95, IAPWS06 = 6, G12_15 = 1215, G9_12 = 912, &
                                  IF97_METASTABLE_VAPOR = 972, IAPWS95_ICE_IH = 9506, G9_12_ICE_IH = 91206
    integer, parameter, public :: LIQUID = 1, VAPOR = 2, SUPERCRITICAL = 3, HIGH_PRESSURE_ICE = 9, &
                                  LIQUID_VAPOR = 10, ICE_LIQUID = 11, ICE_VAPOR = 12
    integer, parameter, public :: WATER_OK = 0, WATER_INVALID = 1, WATER_OUT_OF_RANGE = 2, &
                                  WATER_COEXISTENCE = 3, WATER_UNSUPPORTED = 4, WATER_NO_CONVERGENCE = 5, &
                                  WATER_METASTABLE = 6, WATER_SINGULAR = 7, WATER_AMBIGUOUS = 8
    real(real64), parameter :: critical_T = 647.096d0, critical_p = 22.064d6, critical_rho = 322d0

    integer, parameter, public :: JAC_T = 1, JAC_P = 2, JAC_RHO = 3, JAC_V = 4, &
                                  JAC_H = 5, JAC_S = 6, JAC_U = 7

    !> Analytic Jacobians in SI units for the available coordinate pairs.
    type, public :: type_water_jacobian
        !> Rows use the JAC_* indices. Columns follow the independent variables in the matrix name.
        real(real64) :: tp(7,2) = 0d0, ph(7,2) = 0d0, ps(7,2) = 0d0
        integer :: status = WATER_INVALID
        !> True when the Jacobian follows a branch at a coexistence endpoint.
        logical :: phase_endpoint = .false.
        !> True for each coordinate pair that uniquely specifies the returned branch.
        logical :: tp_available = .false., ph_available = .false., ps_available = .false.
    end type type_water_jacobian

    type, public :: type_water_state
        type(type_iapws_property) :: properties
        integer :: phase = 0
        integer :: formulation = 0
        integer :: status = WATER_INVALID
        !> Iterations of the single-phase temperature solve, or zero for a direct result.
        integer :: iterations = 0
        real(real64), private :: boundary_dTdp = 0d0, boundary_dvhs_dp(3) = 0d0, delta_vhs(3) = 0d0
        ! Bit i identifies a coexisting phase with enumeration i.
        integer :: coexistence = 0
        real(real64) :: quality = -1d0, solid_fraction = -1d0
    contains
        procedure, public :: jacobian => state_jacobian
    end type type_water_state

    type, public :: type_water_saturation
        type(type_iapws_property) :: liquid, vapor
        real(real64) :: T = 0d0, p = 0d0, sigma = 0d0
        integer :: status = WATER_INVALID
        integer :: formulation = 0
    end type type_water_saturation

    type, public :: type_water
        private
        integer :: fluid = 0
        logical :: extend_cold_liquid = .false.
        type(type_iapws95) :: scientific
        type(type_iapws97) :: industrial
        type(type_iapws06) :: ice
    contains
        procedure, public :: state_tp
        procedure, public :: state_tp_checked
        procedure, public :: state_ph_checked
        procedure, public :: state_ps_checked
        procedure, public :: saturation_t
        procedure, public :: saturation_p
        procedure, public :: state_ph
        procedure, public :: state_ps
        procedure, public :: transport
    end type type_water

contains

    pure function new_water(fluid, extend_cold_liquid) result(water)
        implicit none
        integer, intent(in) :: fluid
        logical, intent(in), optional :: extend_cold_liquid
        type(type_water) :: water
        water%fluid = fluid
        if (present(extend_cold_liquid)) water%extend_cold_liquid = extend_cold_liquid
        select case (fluid)
        case (IF97)
            call water%industrial%initialize()
            if (water%extend_cold_liquid) call water%scientific%initialize()
        case (IAPWS95)
            call water%scientific%initialize()
        end select
        call water%ice%initialize()
    end function new_water

    pure elemental subroutine saturation_t(self, T, sat)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: T
        type(type_water_saturation), intent(out) :: sat
        real(real64) :: rl, rv
        logical :: valid

        sat%T = T
        sat%p = ieee_value(0d0, ieee_quiet_nan)
        sat%sigma = sat%p
        sat%liquid = missing_properties(T, sat%p)
        sat%vapor = sat%liquid
        if (.not. ieee_is_finite(T)) return
        sat%status = WATER_OUT_OF_RANGE
        if (T < triple_T .or. T > critical_T) return
        if (T == critical_T) then
            if (self%fluid /= IF97 .and. self%fluid /= IAPWS95) then
                sat%status = WATER_INVALID
                return
            end if
            sat%p = critical_p
            sat%formulation = self%fluid
            sat%sigma = 0d0
            sat%liquid = missing_properties(T, critical_p)
            sat%liquid%rho = critical_rho
            sat%liquid%nu = 1d0 / critical_rho
            sat%vapor = sat%liquid
            sat%status = WATER_SINGULAR
            return
        end if
        select case (self%fluid)
        case (IF97)
            call self%industrial%saturation_endpoints(T, sat%liquid, sat%vapor, valid)
            if (valid) sat%p = sat%liquid%p
        case (IAPWS95)
            call self%scientific%saturation_densities(T, sat%p, rl, rv, sat%liquid, sat%vapor)
            valid = ieee_is_finite(sat%p)
        case default
            sat%status = WATER_INVALID
            return
        end select
        if (valid) then
            sat%formulation = self%fluid
            sat%sigma = surface_tension(T)
            sat%status = WATER_OK
        else
            sat%liquid = missing_properties(T, sat%p)
            sat%vapor = sat%liquid
            sat%status = WATER_NO_CONVERGENCE
        end if
    end subroutine saturation_t

    pure elemental subroutine transport(self, state, mu, conductivity)
        implicit none
        class(type_water), intent(in) :: self
        type(type_water_state), intent(in) :: state
        real(real64), intent(out), optional :: mu, conductivity
        real(real64) :: T, p, rho, nan, m, reference, viscosity_max, conductivity_max
        logical :: need_k, need_mu

        nan = ieee_value(0d0, ieee_quiet_nan)
        if (present(mu)) mu = nan
        if (present(conductivity)) conductivity = nan
        if (state%status /= WATER_OK) return
        if (state%phase < LIQUID .or. state%phase > SUPERCRITICAL) return
        if (state%formulation /= IAPWS95 .and. state%formulation /= IF97) return
        T = state%properties%T
        p = state%properties%p
        rho = state%properties%rho
        if (p < triple_p .and. T < triple_T) return

        viscosity_max = 1173.15d0
        conductivity_max = 1173.15d0
        if (p > 300d6) viscosity_max = 873.15d0
        if (p > 350d6) viscosity_max = 433.15d0
        if (p > 500d6) viscosity_max = 373.15d0
        if (p > 100d6) conductivity_max = 874d0
        if (p > 250d6) conductivity_max = 573d0
        if (p > 687d6) conductivity_max = 403d0
        if (p > 785d6) conductivity_max = 348d0
        need_k = present(conductivity) .and. T <= conductivity_max
        need_mu = present(mu) .and. T <= viscosity_max
        if (.not. need_k .and. .not. need_mu) return

        if (state%formulation == IAPWS95) then
            if (.not. self%scientific%is_initialized) return
            call self%scientific%calc_drho_dp(1.5d0 * critical_T, rho, reference)
            m = viscosity(T, rho, state%properties%drho_dP, reference)
            if (need_k) conductivity = thermal_conductivity(T, rho, state%properties%cp, state%properties%cv, &
                                                           m, state%properties%drho_dP, reference)
        else
            m = viscosity(T, rho)
            if (need_k) conductivity = thermal_conductivity(T, rho, state%properties%cp, state%properties%cv, &
                                                           m, state%properties%drho_dP)
        end if
        if (need_mu) mu = m
    end subroutine transport

    !> Evaluate the analytic derivatives of the returned phase or equilibrium mixture.
    pure elemental subroutine state_jacobian(self, jac)
        implicit none
        class(type_water_state), intent(in) :: self
        type(type_water_jacobian), intent(out) :: jac
        type(type_iapws_property) :: properties
        real(real64) :: dv_dT, dv_dp, dh_dp, ds_dp

        jac%status = self%status
        if (self%status /= WATER_OK .and. self%status /= WATER_METASTABLE) return
        jac%status = WATER_UNSUPPORTED
        select case (self%phase)
        case (LIQUID_VAPOR, ICE_LIQUID, ICE_VAPOR)
            call mixture_jacobian(self, jac)
            return
        case (LIQUID, VAPOR, SUPERCRITICAL, ICE_IH)
            continue
        case default
            return
        end select
        properties = self%properties
        jac%status = WATER_SINGULAR
        if (.not. all(ieee_is_finite([properties%T, properties%p, properties%rho, properties%nu, properties%cp, &
                                     properties%drho_dT, properties%drho_dP]))) return
        if (properties%T <= 0d0 .or. properties%cp <= 0d0 .or. properties%rho <= 0d0) return
        dv_dT = -properties%drho_dT * properties%nu**2
        dv_dp = -properties%drho_dP * properties%nu**2
        dh_dp = properties%nu - properties%T * dv_dT
        ds_dp = -dv_dT
        jac%tp(:,1) = [1d0, 0d0, properties%drho_dT, dv_dT, properties%cp, &
                       properties%cp / properties%T, properties%cp - properties%p * dv_dT]
        jac%tp(:,2) = [0d0, 1d0, properties%drho_dP, dv_dp, dh_dp, ds_dp, -properties%T * dv_dT - properties%p * dv_dp]
        jac%ph(:,2) = jac%tp(:,1) / properties%cp
        jac%ph(:,1) = jac%tp(:,2) - dh_dp * jac%ph(:,2)
        jac%ps(:,2) = properties%T * jac%ph(:,2)
        jac%ps(:,1) = jac%tp(:,2) - ds_dp * jac%ps(:,2)
        if (.not. all(ieee_is_finite(jac%tp)) .or. .not. all(ieee_is_finite(jac%ph)) .or. &
            .not. all(ieee_is_finite(jac%ps))) then
            jac = type_water_jacobian(status=WATER_SINGULAR)
            return
        end if
        jac%tp_available = .true.
        jac%ph_available = .true.
        jac%ps_available = .true.
        jac%phase_endpoint = self%coexistence /= 0
        jac%status = self%status
    end subroutine state_jacobian

    pure subroutine mixture_jacobian(self, jac)
        implicit none
        class(type_water_state), intent(in) :: self
        type(type_water_jacobian), intent(inout) :: jac
        real(real64) :: pressure_column(7), fraction_column(7), dv_dp, dh_dp, ds_dp, delta_v, delta_h, delta_s

        jac%status = WATER_SINGULAR
        associate(properties => self%properties)
            dv_dp = self%boundary_dvhs_dp(1)
            dh_dp = self%boundary_dvhs_dp(2)
            ds_dp = self%boundary_dvhs_dp(3)
            delta_v = self%delta_vhs(1)
            delta_h = self%delta_vhs(2)
            delta_s = self%delta_vhs(3)
            if (delta_h == 0d0 .or. delta_s == 0d0) return
            pressure_column = [self%boundary_dTdp, 1d0, -properties%rho**2 * dv_dp, dv_dp, dh_dp, ds_dp, &
                               dh_dp - properties%nu - properties%p * dv_dp]
            fraction_column = [0d0, 0d0, -properties%rho**2 * delta_v, delta_v, delta_h, delta_s, delta_h - properties%p * delta_v]
        end associate
        jac%ph(:,2) = fraction_column / delta_h
        jac%ph(:,1) = pressure_column - dh_dp * jac%ph(:,2)
        jac%ps(:,2) = fraction_column / delta_s
        jac%ps(:,1) = pressure_column - ds_dp * jac%ps(:,2)
        if (.not. all(ieee_is_finite(jac%ph)) .or. .not. all(ieee_is_finite(jac%ps))) then
            jac = type_water_jacobian(status=WATER_SINGULAR)
            return
        end if
        jac%ph_available = .true.
        jac%ps_available = .true.
        jac%phase_endpoint = self%quality == 0d0 .or. self%quality == 1d0 .or. &
                             self%solid_fraction == 0d0 .or. self%solid_fraction == 1d0
        jac%status = WATER_OK
    end subroutine mixture_jacobian

    pure elemental function bounded_temperature_guess(default_T, lo, hi, guess_T) result(T)
        implicit none
        real(real64) :: T
        real(real64), intent(in) :: default_T, lo, hi
        real(real64), intent(in), optional :: guess_T
        T = default_T
        if (.not. present(guess_T)) return
        if (.not. ieee_is_finite(guess_T)) return
        if (guess_T > lo .and. guess_T < hi) T = guess_T
    end function bounded_temperature_guess

    !> Evaluate a temperature and pressure state with finite failure outputs.
    pure elemental subroutine state_tp_checked(self, T, p, state, phase, allow_metastable, guess_rho, jac)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: T, p
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: phase
        logical, intent(in), optional :: allow_metastable
        real(real64), intent(in), optional :: guess_rho
        type(type_water_jacobian), intent(out), optional :: jac
        call self%state_tp(T, p, state, phase, allow_metastable, guess_rho)
        call finite_solver_result(state, jac)
    end subroutine state_tp_checked

    !> Invert pressure and enthalpy with an optional temperature guess and Jacobian.
    pure elemental subroutine state_ph_checked(self, p, h, state, region, guess_T, jac)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, h
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: region
        real(real64), intent(in), optional :: guess_T
        type(type_water_jacobian), intent(out), optional :: jac
        call self%state_ph(p, h, state, region, guess_T)
        call finite_solver_result(state, jac)
    end subroutine state_ph_checked

    !> Invert pressure and entropy with an optional temperature guess and Jacobian.
    pure elemental subroutine state_ps_checked(self, p, s, state, region, guess_T, jac)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, s
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: region
        real(real64), intent(in), optional :: guess_T
        type(type_water_jacobian), intent(out), optional :: jac
        call self%state_ps(p, s, state, region, guess_T)
        call finite_solver_result(state, jac)
    end subroutine state_ps_checked

    pure elemental subroutine finite_solver_result(state, jacobian)
        implicit none
        type(type_water_state), intent(inout) :: state
        type(type_water_jacobian), intent(out), optional :: jacobian
        type(type_water_jacobian) :: jac
        call state%jacobian(jac)
        if (present(jacobian)) jacobian = jac
        ! Keep equilibrium properties on mixtures and clear undefined single-phase responses.
        if (jac%status == WATER_OK .or. jac%status == WATER_METASTABLE) then
            associate(q => state%properties)
                if (.not. jac%tp_available) then
                    if (all(ieee_is_finite([q%T, q%p, q%nu, q%rho, q%h, q%s, q%u]))) then
                        q = type_iapws_property(region_id=q%region_id, T=q%T, p=q%p, nu=q%nu, rho=q%rho, &
                                                h=q%h, s=q%s, u=q%u)
                        return
                    end if
                end if
                if (all(ieee_is_finite([q%h, q%s, q%u, q%w, q%cv, q%alpha, q%beta, q%kappa_T, q%kappa_s]))) return
            end associate
            jac%status = WATER_SINGULAR
        end if
        state%properties = type_iapws_property()
        state%status = jac%status
        if (present(jacobian)) jacobian = type_water_jacobian(status=state%status)
    end subroutine finite_solver_result

    pure elemental subroutine state_ph(self, p, h, state, region, guess_T)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, h
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: region
        real(real64), intent(in), optional :: guess_T
        logical :: accepted
        if (present(guess_T)) then
            call try_warm_inverse(self, p, h, .false., state, guess_T, region, accepted)
            if (accepted) return
        end if
        if (self%fluid == IAPWS95) then
            if (present(region)) then
                state%properties = missing_properties(ieee_value(0d0, ieee_quiet_nan), p)
                state%status = WATER_INVALID
                return
            end if
            call invert_equilibrium(self, p, h, .false., state, guess_T)
            if (present(guess_T) .and. state%status == WATER_NO_CONVERGENCE) &
                call invert_equilibrium(self, p, h, .false., state)
        else
            call invert_if97(self, p, h, .false., state, region, guess_T)
            if (present(guess_T) .and. state%status == WATER_NO_CONVERGENCE) &
                call invert_if97(self, p, h, .false., state, region)
        end if
    end subroutine state_ph

    pure elemental subroutine state_ps(self, p, s, state, region, guess_T)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, s
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: region
        real(real64), intent(in), optional :: guess_T
        logical :: accepted
        if (present(guess_T)) then
            call try_warm_inverse(self, p, s, .true., state, guess_T, region, accepted)
            if (accepted) return
        end if
        if (self%fluid == IAPWS95) then
            if (present(region)) then
                state%properties = missing_properties(ieee_value(0d0, ieee_quiet_nan), p)
                state%status = WATER_INVALID
                return
            end if
            call invert_equilibrium(self, p, s, .true., state, guess_T)
            if (present(guess_T) .and. state%status == WATER_NO_CONVERGENCE) &
                call invert_equilibrium(self, p, s, .true., state)
        else
            call invert_if97(self, p, s, .true., state, region, guess_T)
            if (present(guess_T) .and. state%status == WATER_NO_CONVERGENCE) &
                call invert_if97(self, p, s, .true., state, region)
        end if
    end subroutine state_ps

    pure subroutine try_warm_inverse(self, p, target, entropy, state, guess_T, region, accepted)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, target, guess_T
        logical, intent(in) :: entropy
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: region
        logical, intent(out) :: accepted
        type(type_water_state) :: candidate
        real(real64) :: T, target_dT, temperature_step, residual
        integer :: iteration, initial_phase

        accepted = .false.
        if (.not. all(ieee_is_finite([p, target, guess_T]))) return
        if (p < 1d-100 .or. p > 1d9 .or. abs(target) > 1d9) return
        if (self%fluid /= IF97 .and. self%fluid /= IAPWS95) return
        if (present(region)) then
            if (self%fluid /= IF97 .or. region < 1 .or. region > 5 .or. region == 4) return
        end if
        T = guess_T
        initial_phase = 0
        do iteration = 1, 4
            if (self%fluid == IF97 .and. .not. present(region)) then
                ! Leave inversion across potentially overlapping IF97 joins to the full solver.
                if (p > 16d6 .or. T >= 1070d0) return
            end if
            call self%state_tp(T, p, candidate)
            if (candidate%status /= WATER_OK .or. candidate%coexistence /= 0) return
            if (self%fluid == IF97 .and. candidate%formulation /= IF97) return
            if (present(region)) then
                if (candidate%properties%region_id /= region) return
            end if
            if (iteration == 1) initial_phase = candidate%phase
            if (candidate%phase /= initial_phase) return
            target_dT = candidate%properties%cp
            if (.not. ieee_is_finite(target_dT)) return
            if (target_dT <= 1d0) return
            if (entropy) target_dT = target_dT / T
            residual = merge(candidate%properties%s, candidate%properties%h, entropy) - target
            if (abs(residual) <= 32d0 * epsilon(T) * max(abs(target), abs(target_dT * T))) then
                state = candidate
                state%iterations = iteration
                accepted = .true.
                return
            end if
            temperature_step = residual / target_dT
            if (abs(temperature_step) > 0.025d0 * T) return
            T = T - temperature_step
        end do
    end subroutine try_warm_inverse

    pure subroutine invert_equilibrium(self, p, target, entropy, state, guess_T)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, target
        logical, intent(in) :: entropy
        real(real64), intent(in), optional :: guess_T
        type(type_water_state), intent(out) :: state
        type(type_water_state) :: endpoint
        type(type_water_saturation) :: sat
        type(type_iapws_property) :: lower, upper, solid, fluid
        real(real64) :: lo, hi, Tm, Ts, rho, a, b, mid, residual, x
        integer :: i, phase_above, melting_solid
        logical :: has_ice

        state%properties = missing_properties(ieee_value(0d0, ieee_quiet_nan), p)
        if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(target)) return
        if (p < 1d-100) return
        state%status = WATER_OUT_OF_RANGE
        if (p > 1d9) return
        lo = 50d0
        hi = 1273d0
        has_ice = p >= sublimation_pressure(lo)
        if (has_ice) then
            a = lo
            b = triple_T
            phase_above = VAPOR
            if (p >= triple_p) then
                a = 251.165d0
                phase_above = LIQUID
            end if
            melting_solid = ICE_IH
            if (p > 632.4d6) then
                melting_solid = ICE_VI
                a = 273.31d0
                b = 355d0
            else if (p > 350.1d6) then
                melting_solid = ICE_V
                a = 256.164d0
                b = 273.31d0
            else if (p > 208.566d6) then
                melting_solid = ICE_III
                a = 251.165d0
                b = 256.164d0
            end if
            do i = 1, 54
                mid = (a + b) * 0.5d0
                if (phase_above == VAPOR) then
                    residual = sublimation_pressure(mid) - p
                else
                    residual = melting_pressure(mid, melting_solid) - p
                    if (melting_solid == ICE_IH) residual = -residual
                end if
                if (residual < 0d0) then
                    a = mid
                else
                    b = mid
                end if
            end do
            Tm = (a + b) * 0.5d0
            if (p == triple_p) then
                call self%saturation_t(triple_T, sat)
                rho = sat%liquid%rho
            else
                call self%scientific%calc_rho(Tm, p, rho)
            end if
            state%status = WATER_NO_CONVERGENCE
            if (.not. ieee_is_finite(rho)) return
            call self%scientific%calc_properties(Tm, rho, fluid)
            fluid%p = p
            b = merge(fluid%s, fluid%h, entropy)
            if (target < b .and. melting_solid /= ICE_IH) then
                state%phase = HIGH_PRESSURE_ICE
                state%status = WATER_UNSUPPORTED
                return
            end if
            if (target <= b .and. melting_solid == ICE_IH) then
                call self%ice%calc_properties(Tm, p, solid)
                a = merge(solid%s, solid%h, entropy)
                if (target >= a) then
                    if (p == triple_p) then
                        state%status = WATER_COEXISTENCE
                        state%properties = missing_properties(Tm, p)
                        state%coexistence = ibset(ibset(ibset(0, ICE_IH), LIQUID), VAPOR)
                        return
                    end if
                    x = (target - a) / (b - a)
                    if (phase_above == LIQUID) then
                        call set_mixture(solid, fluid, x, p, ICE_IH, LIQUID, ICE_LIQUID, IAPWS95_ICE_IH, state, &
                                         1d0 / ice_ih_melting_slope(Tm))
                    else
                        call set_mixture(solid, fluid, x, p, ICE_IH, VAPOR, ICE_VAPOR, &
                                         merge(G9_12_ICE_IH, IAPWS95_ICE_IH, Tm < 130d0), state, 1d0 / sublimation_slope(Tm))
                    end if
                    state%solid_fraction = 1d0 - x
                    return
                end if
                hi = Tm
                upper = solid
                call self%ice%calc_properties(lo, p, lower)
                call solve_segment(self, p, target, entropy, lo, hi, lower, upper, state, guess_T)
                return
            end if
            lo = Tm
            lower = fluid
        else
            call self%state_tp(lo, p, endpoint)
            if (endpoint%status /= WATER_OK) return
            lower = endpoint%properties
        end if

        call self%state_tp(hi, p, endpoint)
        if (endpoint%status /= WATER_OK) return
        upper = endpoint%properties
        if (p > triple_p .and. p < 22.064d6) then
            call saturation_at_pressure95(self, p, sat, Ts)
            state%status = sat%status
            if (sat%status == WATER_OK) then
                a = merge(sat%liquid%s, sat%liquid%h, entropy)
                b = merge(sat%vapor%s, sat%vapor%h, entropy)
                if (target >= a .and. target <= b) then
                    x = (target - a) / (b - a)
                    call set_mixture(sat%liquid, sat%vapor, x, p, LIQUID, VAPOR, LIQUID_VAPOR, IAPWS95, state, &
                        Ts * (sat%vapor%nu - sat%liquid%nu) / (sat%vapor%h - sat%liquid%h))
                    state%quality = x
                    return
                else if (target < a) then
                    hi = Ts
                    upper = sat%liquid
                else
                    lo = Ts
                    lower = sat%vapor
                end if
            else
                call self%state_tp(max(lo, Ts - 0.125d0), p, endpoint)
                if (endpoint%status == WATER_OK) then
                    if (target <= merge(endpoint%properties%s, endpoint%properties%h, entropy)) then
                        hi = endpoint%properties%T
                        upper = endpoint%properties
                        call solve_segment(self, p, target, entropy, lo, hi, lower, upper, state, guess_T)
                        return
                    end if
                end if
                call self%state_tp(min(hi, Ts + 0.125d0), p, endpoint)
                if (endpoint%status == WATER_OK) then
                    if (target >= merge(endpoint%properties%s, endpoint%properties%h, entropy)) then
                        lo = endpoint%properties%T
                        lower = endpoint%properties
                        call solve_segment(self, p, target, entropy, lo, hi, lower, upper, state, guess_T)
                    end if
                end if
                return
            end if
        end if
        if (p == triple_p) then
            if (target <= merge(sat%vapor%s, sat%vapor%h, entropy)) then
                state%status = WATER_COEXISTENCE
                state%properties = missing_properties(triple_T, p)
                state%coexistence = ibset(ibset(ibset(0, ICE_IH), LIQUID), VAPOR)
                return
            end if
            lo = triple_T
            lower = sat%vapor
        end if
        call solve_segment(self, p, target, entropy, lo, hi, lower, upper, state, guess_T)
    end subroutine invert_equilibrium

    pure subroutine saturation_at_pressure95(self, p, sat, Ts)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p
        type(type_water_saturation), intent(out) :: sat
        real(real64), intent(out) :: Ts
        real(real64) :: a, b, residual, slope, trial
        integer :: i

        a = triple_T
        b = critical_T
        Ts = (a + b) * 0.5d0
        do i = 1, 22
            if (self%scientific%saturation_estimate(Ts) < p) then
                a = Ts
            else
                b = Ts
            end if
            Ts = (a + b) * 0.5d0
        end do
        a = triple_T
        b = critical_T
        do i = 1, 16
            call self%saturation_t(Ts, sat)
            if (sat%status /= WATER_OK) return
            residual = sat%p - p
            if (abs(residual) <= max(5d-9, 128d0 * epsilon(p) * p)) return
            if (residual > 0d0) then
                b = Ts
            else
                a = Ts
            end if
            slope = (sat%vapor%h - sat%liquid%h) / (Ts * (sat%vapor%nu - sat%liquid%nu))
            trial = Ts - residual / slope
            if (trial <= a .or. trial >= b) trial = (a + b) * 0.5d0
            if (trial == Ts) exit
            Ts = trial
        end do
        sat%status = WATER_NO_CONVERGENCE
    end subroutine saturation_at_pressure95

    pure subroutine set_mixture(lower, upper, upper_fraction, p, lower_phase, upper_phase, mixture, formulation, state, dTdp)
        implicit none
        type(type_iapws_property), intent(in) :: lower, upper
        real(real64), intent(in) :: upper_fraction, p, dTdp
        integer, intent(in) :: lower_phase, upper_phase, mixture, formulation
        type(type_water_state), intent(out) :: state
        state%status = WATER_OK
        state%phase = mixture
        state%formulation = formulation
        state%coexistence = ibset(ibset(0, lower_phase), upper_phase)
        state%properties = missing_properties(lower%T, p)
        state%properties%nu = (1d0 - upper_fraction) * lower%nu + upper_fraction * upper%nu
        state%properties%rho = 1d0 / state%properties%nu
        state%properties%h = (1d0 - upper_fraction) * lower%h + upper_fraction * upper%h
        state%properties%s = (1d0 - upper_fraction) * lower%s + upper_fraction * upper%s
        state%properties%u = state%properties%h - p * state%properties%nu
        state%boundary_dTdp = dTdp
        state%boundary_dvhs_dp = (1d0 - upper_fraction) * endpoint_boundary_derivatives(lower, dTdp) + &
                                          upper_fraction * endpoint_boundary_derivatives(upper, dTdp)
        state%delta_vhs = [upper%nu - lower%nu, upper%h - lower%h, upper%s - lower%s]
    end subroutine set_mixture

    pure function endpoint_boundary_derivatives(properties, dTdp) result(dvhs)
        implicit none
        type(type_iapws_property), intent(in) :: properties
        real(real64), intent(in) :: dTdp
        real(real64) :: dvhs(3), dv_dT
        dv_dT = -properties%drho_dT * properties%nu**2
        dvhs = [-properties%drho_dP * properties%nu**2 + dv_dT * dTdp, &
                properties%nu - properties%T * dv_dT + properties%cp * dTdp, -dv_dT + (properties%cp / properties%T) * dTdp]
    end function endpoint_boundary_derivatives

    pure subroutine solve_segment(self, p, target, entropy, low_T, high_T, lower, upper, state, guess_T)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, target, low_T, high_T
        logical, intent(in) :: entropy
        type(type_iapws_property), intent(in) :: lower, upper
        real(real64), intent(in), optional :: guess_T
        type(type_water_state), intent(inout) :: state
        type(type_water_state) :: candidate
        real(real64) :: lo, hi, a, b, T, trial, residual, slope
        integer :: i

        state%status = WATER_OUT_OF_RANGE
        a = merge(lower%s, lower%h, entropy)
        b = merge(upper%s, upper%h, entropy)
        if (target < a .or. target > b .or. b <= a) return
        lo = low_T
        hi = high_T
        T = lo + (hi - lo) * ((target - a) / (b - a))
        T = bounded_temperature_guess(T, lo, hi, guess_T)
        state%status = WATER_NO_CONVERGENCE
        do i = 1, 80
            call self%state_tp(T, p, candidate)
            if (candidate%status /= WATER_OK) return
            residual = merge(candidate%properties%s, candidate%properties%h, entropy) - target
            slope = candidate%properties%cp
            if (entropy) slope = slope / T
            if (abs(residual) <= 256d0 * epsilon(T) * max(abs(target), abs(slope * T))) then
                state = candidate
                state%iterations = i
                return
            end if
            if (residual > 0d0) then
                hi = T
            else
                lo = T
            end if
            trial = lo + (hi - lo) * 0.5d0
            if (slope > 0d0) then
                if (T - residual / slope > lo .and. T - residual / slope < hi) trial = T - residual / slope
            end if
            if (trial == T) return
            T = trial
        end do
    end subroutine solve_segment

    pure subroutine invert_if97(self, p, target, entropy, state, region_hint, guess_T)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p, target
        logical, intent(in) :: entropy
        real(real64), intent(in), optional :: guess_T
        type(type_water_state), intent(out) :: state
        type(type_water_saturation) :: sat
        integer, intent(in), optional :: region_hint
        real(real64) :: joins(3), edge_a, edge_b
        type(type_iapws_property) :: lower, upper, properties
        real(real64) :: lo, hi, Ts, T, a, b, residual, slope, trial, tolerance, x
        integer :: iteration, region

        state%properties = missing_properties(ieee_value(0d0, ieee_quiet_nan), p)
        if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(target)) return
        if (p < 1d-100) return
        if (present(region_hint)) then
            if (region_hint < 1 .or. region_hint > 5) return
        end if
        state%status = WATER_UNSUPPORTED
        if (self%fluid /= IF97) return
        state%status = WATER_OUT_OF_RANGE
        if (p > 100d6) return
        lo = 273.15d0
        hi = merge(1073.15d0, 2273.15d0, p > 50d6)
        Ts = self%industrial%saturation_temperature(p)
        call self%industrial%calc_properties(lo, p, lower)
        call self%industrial%calc_properties(hi, p, upper)
        if (ieee_is_finite(Ts)) then
            if (Ts >= triple_T .and. Ts < critical_T) then
                call self%saturation_t(Ts, sat)
                state%status = sat%status
                if (sat%status /= WATER_OK) then
                    trial = max(lo, Ts - 0.01d0)
                    call self%industrial%calc_properties(trial, p, properties)
                    if (target <= merge(properties%s, properties%h, entropy)) then
                        hi = trial
                        upper = properties
                    else
                        trial = min(hi, Ts + 0.01d0)
                        call self%industrial%calc_properties(trial, p, properties)
                        if (target < merge(properties%s, properties%h, entropy)) return
                        lo = trial
                        lower = properties
                    end if
                else
                    a = merge(sat%liquid%s, sat%liquid%h, entropy)
                    b = merge(sat%vapor%s, sat%vapor%h, entropy)
                    if (target >= a .and. target <= b) then
                        if (present(region_hint)) then
                            if (region_hint /= 4) then
                                state%status = WATER_OUT_OF_RANGE
                                return
                            end if
                        end if
                        x = (target - a) / (b - a)
                        call set_mixture(sat%liquid, sat%vapor, x, p, LIQUID, VAPOR, LIQUID_VAPOR, IF97, state, &
                                         1d0 / self%industrial%saturation_slope(Ts))
                        state%quality = x
                        state%properties%region_id = 4
                        return
                    else if (target < a) then
                        hi = Ts
                        upper = sat%liquid
                    else
                        lo = Ts
                        lower = sat%vapor
                    end if
                end if
            end if
        end if
        joins = [623.15d0, 1073.15d0, 0d0]
        if (p > 16.5291642526d6) then
            joins(3) = 572.54459862746d0 + sqrt((p / 1d6 - 13.918839778870d0) / 0.0010192970039326d0)
        end if
        if (present(region_hint)) then
            trial = lo
            T = hi
            select case (region_hint)
            case (1)
                T = min(hi, joins(1))
            case (2)
                if (joins(3) > 0d0) trial = max(lo, joins(3) + 8d0 * spacing(joins(3)))
                T = min(hi, joins(2))
            case (3)
                if (joins(3) <= 0d0) then
                    state%status = WATER_OUT_OF_RANGE
                    return
                end if
                trial = max(lo, nearest(joins(1), 1d0))
                T = min(hi, joins(3) - 8d0 * spacing(joins(3)))
            case (4)
                state%status = WATER_OUT_OF_RANGE
                return
            case (5)
                trial = max(lo, nearest(joins(2), 1d0))
            end select
            if (trial >= T) then
                state%status = WATER_OUT_OF_RANGE
                return
            end if
            if (trial /= lo) call self%industrial%calc_properties(trial, p, lower)
            if (T /= hi) call self%industrial%calc_properties(T, p, upper)
            lo = trial
            hi = T
        else
            do iteration = 1, size(joins)
                trial = joins(iteration)
                if (trial < lo .or. trial > hi) cycle
                if (iteration == 1 .and. joins(3) == 0d0) cycle
                call self%industrial%calc_properties(trial, p, properties)
                edge_a = merge(properties%s, properties%h, entropy)
                if (target == edge_a) then
                    call self%state_tp(trial, p, state)
                    return
                end if
                if (trial - lo <= 16d0 * spacing(trial) .or. hi - trial <= 16d0 * spacing(trial)) cycle
                if (abs(target - edge_a) > 0.002d0 * max(1d0, abs(edge_a))) cycle
                call self%industrial%calc_properties(trial - 8d0 * spacing(trial), p, properties)
                edge_a = merge(properties%s, properties%h, entropy)
                call self%industrial%calc_properties(trial + 8d0 * spacing(trial), p, properties)
                edge_b = merge(properties%s, properties%h, entropy)
                if (edge_b < edge_a .and. target > edge_b .and. target < edge_a) then
                    state%status = WATER_AMBIGUOUS
                    return
                end if
            end do
        end if
        state%status = WATER_OUT_OF_RANGE
        a = merge(lower%s, lower%h, entropy)
        b = merge(upper%s, upper%h, entropy)
        if (target < a .or. target > b) return
        T = lo + (hi - lo) * ((target - a) / (b - a))
        region = self%industrial%get_region(T, p)
        if (p > 16.5291642526d6) then
            trial = 623.15d0
            call self%industrial%calc_properties(trial, p, properties)
            if (target <= merge(properties%s, properties%h, entropy)) then
                region = 1
            else
                trial = 572.54459862746d0 + sqrt((p * 1d-6 - 13.918839778870d0) / 0.0010192970039326d0)
                call self%industrial%calc_properties(trial, p, properties)
                region = merge(3, 2, target < merge(properties%s, properties%h, entropy))
            end if
        else if (ieee_is_finite(Ts)) then
            region = merge(1, 2, hi == Ts)
        end if
        if (present(region_hint)) region = region_hint
        trial = self%industrial%backward_temperature(p, target, region, entropy)
        if (ieee_is_finite(trial)) then
            if (trial > lo .and. trial < hi) T = trial
        end if
        T = bounded_temperature_guess(T, lo, hi, guess_T)
        state%status = WATER_NO_CONVERGENCE
        do iteration = 1, 80
            call self%industrial%calc_properties(T, p, properties)
            if (.not. ieee_is_finite(properties%rho)) return
            residual = merge(properties%s, properties%h, entropy) - target
            slope = properties%cp
            if (entropy) slope = slope / T
            tolerance = 32d0 * epsilon(T) * max(abs(target), abs(slope * T))
            if (abs(residual) <= tolerance) then
                call self%state_tp(T, p, state)
                state%iterations = iteration
                if (present(region_hint)) then
                    if (state%properties%region_id /= region_hint) then
                        state%properties = missing_properties(T, p)
                        state%formulation = 0
                        state%status = WATER_OUT_OF_RANGE
                        return
                    end if
                end if
                if (state%phase == ICE_IH .or. state%phase == HIGH_PRESSURE_ICE) then
                    state%properties = missing_properties(T, p)
                    state%formulation = 0
                    state%status = WATER_UNSUPPORTED
                end if
                return
            end if
            if (residual > 0d0) then
                hi = T
            else
                lo = T
            end if
            trial = lo + (hi - lo) * 0.5d0
            if (slope > 0d0) then
                if (T - residual / slope > lo .and. T - residual / slope < hi) trial = T - residual / slope
            end if
            if (trial == T) return
            T = trial
        end do
    end subroutine invert_if97

    pure elemental subroutine saturation_p(self, p, sat)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: p
        type(type_water_saturation), intent(out) :: sat
        real(real64) :: T, nan

        nan = ieee_value(0d0, ieee_quiet_nan)
        sat%T = nan
        sat%p = p
        sat%sigma = nan
        sat%liquid = missing_properties(nan, p)
        sat%vapor = sat%liquid
        if (.not. ieee_is_finite(p)) return
        if (p < 1d-100) return
        sat%status = WATER_OUT_OF_RANGE
        if (p > critical_p) return
        if (p == critical_p) then
            call self%saturation_t(critical_T, sat)
            return
        end if
        select case (self%fluid)
        case (IF97)
            T = self%industrial%saturation_temperature(p)
            if (.not. ieee_is_finite(T)) return
            if (T < triple_T) then
                if (on_boundary(p, self%industrial%saturation_pressure(triple_T))) then
                    T = triple_T
                else
                    return
                end if
            end if
            call self%saturation_t(T, sat)
        case (IAPWS95)
            if (p < 611.65d0) return
            if (p < 611.66d0) then
                call self%saturation_t(triple_T, sat)
                if (sat%status /= WATER_OK) return
                if (on_boundary(p, sat%p)) return
                if (p < sat%p) then
                    sat%status = WATER_OUT_OF_RANGE
                    sat%T = nan
                    sat%p = p
                    sat%sigma = nan
                    sat%formulation = 0
                    sat%liquid = missing_properties(nan, p)
                    sat%vapor = sat%liquid
                    return
                end if
            end if
            call saturation_at_pressure95(self, p, sat, T)
        case default
            sat%status = WATER_INVALID
        end select
    end subroutine saturation_p

    pure elemental function missing_properties(T, p) result(properties)
        implicit none
        real(real64), intent(in) :: T, p
        type(type_iapws_property) :: properties
        real(real64) :: nan
        nan = ieee_value(0d0, ieee_quiet_nan)
        properties = type_iapws_property(-1, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, &
                                         p, T, nan, nan, nan, nan)
    end function missing_properties

    pure elemental subroutine state_tp(self, T, p, state, phase, allow_metastable, guess_rho)
        implicit none
        class(type_water), intent(in) :: self
        real(real64), intent(in) :: T, p
        type(type_water_state), intent(out) :: state
        integer, intent(in), optional :: phase
        logical, intent(in), optional :: allow_metastable
        real(real64), intent(in), optional :: guess_rho
        real(real64) :: boundary, rho, rho_liq, rho_vap
        type(type_water_saturation) :: sat
        logical :: scientific, valid, solid_stable, metastable

        state%properties = missing_properties(T, p)
        if (.not. ieee_is_finite(T) .or. .not. ieee_is_finite(p)) return
        if (T <= 0d0 .or. p < 1d-100) return
        if (self%fluid /= IF97 .and. self%fluid /= IAPWS95) return
        state%status = WATER_OUT_OF_RANGE
        if (T < 50d0 .or. p > 1d9) return
        if (T == critical_T .and. p == critical_p) then
            state%status = WATER_SINGULAR
            state%phase = SUPERCRITICAL
            state%formulation = self%fluid
            state%properties%rho = critical_rho
            state%properties%nu = 1d0 / critical_rho
            return
        end if

        solid_stable = .false.
        if (T <= triple_T .and. p <= 208.566d6) then
            state%phase = ICE_IH
            if (p <= triple_p) then
                boundary = sublimation_pressure(T)
                if (on_boundary(p, boundary)) state%coexistence = ibset(ibset(0, ICE_IH), VAPOR)
                if (p < boundary) state%phase = VAPOR
            else if (T >= 251.165d0) then
                boundary = melting_pressure(T)
                if (on_boundary(p, boundary)) state%coexistence = ibset(ibset(0, ICE_IH), LIQUID)
                if (p > boundary) state%phase = LIQUID
            end if
            if (T == triple_T .and. on_boundary(p, triple_p)) then
                state%coexistence = ibset(ibset(ibset(0, ICE_IH), VAPOR), LIQUID)
            else if (T == triple_T .and. self%fluid == IAPWS95) then
                if (abs(p - triple_p) < 1d-2) then
                    call self%scientific%saturation_densities(T, boundary, rho_liq, rho_vap)
                    if (ieee_is_finite(boundary)) then
                        if (on_boundary(p, boundary)) state%coexistence = ibset(ibset(0, VAPOR), LIQUID)
                    end if
                end if
            end if
            solid_stable = state%phase == ICE_IH .and. state%coexistence == 0
        else if (p > 208.566d6) then
            boundary = high_pressure_melting_pressure(T)
            if (on_boundary(p, boundary)) then
                state%coexistence = ibset(ibset(0, LIQUID), HIGH_PRESSURE_ICE)
                state%phase = LIQUID
            else if (p > boundary) then
                state%phase = HIGH_PRESSURE_ICE
                solid_stable = .true.
            end if
        end if

        metastable = .false.
        if (present(phase) .and. present(allow_metastable)) then
            metastable = phase == LIQUID .and. allow_metastable .and. solid_stable
        end if
        if (metastable) then
            state%phase = LIQUID
            if (.not. in_supercooled_domain(T, p)) return
            call evaluate_supercooled(T, p, state%properties, valid)
            state%status = WATER_NO_CONVERGENCE
            if (.not. valid) then
                state%properties = missing_properties(T, p)
                return
            end if
            state%formulation = G12_15
            state%status = WATER_METASTABLE
            return
        end if
        if (state%phase == HIGH_PRESSURE_ICE) then
            state%status = WATER_UNSUPPORTED
            return
        end if

        scientific = self%fluid == IAPWS95 .or. &
                     (self%extend_cold_liquid .and. state%phase == LIQUID .and. T < 273.15d0)
        if (state%phase == 0 .and. .not. scientific) then
            if (T >= critical_T) then
                state%phase = merge(SUPERCRITICAL, VAPOR, p >= 22.064d6)
            else
                boundary = self%industrial%saturation_pressure(T)
                if (.not. ieee_is_finite(boundary)) return
                if (on_boundary(p, boundary)) state%coexistence = ibset(ibset(0, LIQUID), VAPOR)
                state%phase = merge(LIQUID, VAPOR, p > boundary)
            end if
        end if
        if (state%phase == 0 .and. scientific .and. T < critical_T) then
            boundary = self%scientific%saturation_estimate(T)
            if (abs(p - boundary) < 0.001d0 * boundary) then
                call self%scientific%saturation_densities(T, boundary, rho_liq, rho_vap)
                state%status = WATER_NO_CONVERGENCE
                if (.not. ieee_is_finite(boundary)) return
                if (on_boundary(p, boundary)) state%coexistence = ibset(ibset(0, LIQUID), VAPOR)
            end if
        end if
        if (present(phase) .and. present(allow_metastable)) then
            if (phase == VAPOR .and. allow_metastable .and. self%fluid == IF97 .and. state%phase == LIQUID) then
                call self%industrial%metastable_vapor(T, p, state%properties, valid)
                if (valid) then
                    state%phase = VAPOR
                    state%formulation = IF97_METASTABLE_VAPOR
                    state%status = WATER_METASTABLE
                end if
                return
            end if
        end if
        if (state%coexistence /= 0) then
            if (.not. present(phase)) then
                state%phase = 0
                state%status = WATER_COEXISTENCE
                return
            end if
            if (phase == HIGH_PRESSURE_ICE) then
                if (btest(state%coexistence, HIGH_PRESSURE_ICE)) then
                    state%phase = phase
                    state%status = WATER_UNSUPPORTED
                    return
                end if
            end if
            if (phase < LIQUID .or. phase > ICE_IH) return
            if (.not. btest(state%coexistence, phase)) return
            state%phase = phase
            scientific = self%fluid == IAPWS95 .or. &
                         (self%extend_cold_liquid .and. phase == LIQUID .and. T < 273.15d0)
        else if (present(phase)) then
            if (phase /= state%phase .and. state%phase /= 0) then
                state%status = WATER_UNSUPPORTED
                return
            end if
        end if
        if (state%phase == ICE_IH) then
            call self%ice%calc_properties(T, p, state%properties)
            state%formulation = IAPWS06
        else if (scientific) then
            if (T > 1273d0 .or. p > 1d9) return
            if (btest(state%coexistence, LIQUID) .and. btest(state%coexistence, VAPOR)) then
                if (btest(state%coexistence, ICE_IH)) then
                    call self%scientific%saturation_densities(T, boundary, rho_liq, rho_vap)
                    if (.not. ieee_is_finite(boundary)) then
                        state%status = WATER_NO_CONVERGENCE
                        return
                    end if
                end if
                rho = merge(rho_liq, rho_vap, state%phase == LIQUID)
            else
                call self%scientific%solve_density(T, p, rho, guess_rho)
            end if
            state%status = WATER_NO_CONVERGENCE
            if (.not. ieee_is_finite(rho)) return
            if (state%phase == 0) then
                if (T >= critical_T) then
                    state%phase = merge(SUPERCRITICAL, VAPOR, p >= 22.064d6)
                else
                    state%phase = merge(LIQUID, VAPOR, rho > 322d0)
                end if
                if (present(phase)) then
                    if (phase /= state%phase) then
                        state%status = WATER_UNSUPPORTED
                        return
                    end if
                end if
            end if
            call self%scientific%calc_properties(T, rho, state%properties)
            if (state%coexistence /= 0) state%properties%p = p
            state%formulation = merge(G9_12, IAPWS95, T < 130d0)
        else
            if (T < 273.15d0 .or. self%industrial%get_region(T, p) < 0) return
            if (btest(state%coexistence, LIQUID) .and. btest(state%coexistence, VAPOR)) then
                call self%saturation_t(T, sat)
                state%status = sat%status
                if (sat%status /= WATER_OK) return
                if (state%phase == LIQUID) then
                    state%properties = sat%liquid
                else
                    state%properties = sat%vapor
                end if
                state%properties%p = p
            else
                call self%industrial%evaluate_state(T, p, state%properties, guess_rho)
            end if
            state%formulation = IF97
        end if
        state%status = WATER_NO_CONVERGENCE
        if (.not. ieee_is_finite(state%properties%rho)) return
        state%status = WATER_OK
    end subroutine state_tp

    pure elemental function on_boundary(p, boundary) result(is_on_boundary)
        implicit none
        logical :: is_on_boundary
        real(real64), intent(in) :: p, boundary
        is_on_boundary = abs(p - boundary) <= 32d0 * epsilon(p) * max(abs(boundary), p)
    end function on_boundary

    pure elemental function high_pressure_melting_pressure(T) result(boundary)
        implicit none
        real(real64) :: boundary
        real(real64), intent(in) :: T
        boundary = 0d0
        if (T < 251.165d0) return
        if (T >= 355d0) then
            ! The VII phase lies above the implemented fluid EOS pressure limits.
            boundary = huge(T)
            return
        else if (T >= 273.31d0) then
            boundary = melting_pressure(T, ICE_VI)
        else if (T >= 256.164d0) then
            boundary = melting_pressure(T, ICE_V)
        else
            boundary = melting_pressure(T, ICE_III)
        end if
    end function high_pressure_melting_pressure
end module iapws
