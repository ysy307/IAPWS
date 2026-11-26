submodule(module_iapws97) iapws97_base
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use module_iapws97_backward, only: backward_temperature
    use module_iapws97_volume, only: volume_region3
    implicit none
    !> Invalid region identifier
    integer(int32), parameter :: IAPWS97_INVALID = -1
    !> Region 1 (Liquid water) identifier
    integer(int32), parameter :: IAPWS97_REGION_1 = 1
    !> Region 2 (Steam) identifier
    integer(int32), parameter :: IAPWS97_REGION_2 = 2
    !> Region 3 (Critical / Supercritical) identifier
    integer(int32), parameter :: IAPWS97_REGION_3 = 3
    !> Region 5 (High temperature) identifier
    integer(int32), parameter :: IAPWS97_REGION_5 = 5

    !> Global maximum temperature limit [K]
    real(real64), parameter :: IAPWS97_LIMIT_T_MAX = 2273.15d0
    !> Global maximum pressure limit [Pa]
    real(real64), parameter :: IAPWS97_LIMIT_P_MAX = 100.0d6

    !> Region 5 minimum temperature [K]
    real(real64), parameter :: IAPWS97_R5_T_MIN = 1073.15d0
    !> Region 5 maximum pressure [Pa]
    real(real64), parameter :: IAPWS97_R5_P_MAX = 50.0d6
    ! Phase data shared by saturation values, partials, and total derivatives.
    type :: saturation_phase
        real(real64) :: rho, cp, h, rho_T, rho_p, cp_T, cp_p, h_p
    end type
contains

    module pure elemental subroutine metastable_vapor_iapws97(self, T, p, property, valid)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T, p
        type(type_iapws_property), intent(inout) :: property
        logical, intent(out) :: valid
        type(type_iapws_gibbs_coefficient) :: g
        real(real64) :: Ts, hf, hg, h

        valid = .false.
        if (.not. ieee_is_finite(T) .or. .not. ieee_is_finite(p)) return
        if (T < 273.16d0 .or. p < 611.657d0 .or. p > 10d6) return
        Ts = self%region4%calc_tsat(p)
        if (T > Ts) return
        call self%region2%calc_gamma_metastable(T, p, g)
        h = g%g - T * g%g_t
        call self%region1%calc_h(Ts, p, hf)
        call self%region2%calc_h(Ts, p, hg)
        if (h < hg - 0.05d0 * (hg - hf)) return
        if (g%g_p <= 0d0 .or. g%g_pp >= 0d0 .or. g%g_tt >= 0d0) return
        if (g%g_pp - g%g_pt**2 / g%g_tt >= 0d0) return
        call self%region2%project_properties(T, p, property, g)
        property%region_id = 2
        valid = .true.
    end subroutine metastable_vapor_iapws97

    module pure elemental function backward_volume_iapws97(self, T, p) result(v)
        implicit none
        real(real64) :: v
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T, p
        real(real64) :: Ts
        v = ieee_value(0d0, ieee_quiet_nan)
        if (self%get_region(T, p) /= IAPWS97_REGION_3) return
        Ts = 647.096d0
        if (p <= 22.064d6) Ts = self%region4%calc_tsat(p)
        v = volume_region3(T, p / 1d6, Ts)
    end function backward_volume_iapws97

    module pure elemental function backward_temperature_iapws97(self, p, input, region, entropy) result(T)
        implicit none
        real(real64) :: T
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: p, input
        integer, intent(in) :: region
        logical, intent(in), optional :: entropy
        logical :: use_s
        use_s = .false.
        if (present(entropy)) use_s = entropy
        T = backward_temperature(p, input, region, use_s)
    end function backward_temperature_iapws97

    module pure elemental subroutine saturation_endpoints_iapws97(self, T, liquid, vapor, valid)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T
        type(type_iapws_property), intent(out) :: liquid, vapor
        logical, intent(out) :: valid
        real(real64) :: p, rl, rv

        valid = .false.
        if (.not. ieee_is_finite(T)) return
        if (T < 273.15d0 .or. T >= 647.096d0) return
        p = self%region4%calc_psat(T)
        if (T <= 623.15d0) then
            call self%region1%calc_properties(T, p, liquid)
            call self%region2%calc_properties(T, p, vapor)
            liquid%region_id = 1
            vapor%region_id = 2
        else
            call saturation_root(self, T, p, .false., rl)
            call saturation_root(self, T, p, .true., rv)
            if (.not. ieee_is_finite(rl) .or. .not. ieee_is_finite(rv)) return
            call self%region3%calc_properties(T, rl, liquid)
            call self%region3%calc_properties(T, rv, vapor)
            liquid%region_id = 3
            vapor%region_id = 3
        end if
        liquid%p = p
        vapor%p = p
        valid = .true.
    end subroutine saturation_endpoints_iapws97

    module procedure saturation_slope_iapws97
        implicit none
        slope = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 273.15d0 .or. T > 647.096d0) return
        slope = self%region4%calc_dpsat_dT(T)
    end procedure saturation_slope_iapws97

    module pure elemental function saturation_pressure_iapws97(self, T) result(p_sat)
        implicit none
        real(real64) :: p_sat
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T
        p_sat = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 273.15d0 .or. T > 647.096d0) return
        p_sat = self%region4%calc_psat(T)
        if (T == 647.096d0) p_sat = 22.064d6
    end function saturation_pressure_iapws97

    module pure elemental function saturation_temperature_iapws97(self, p) result(T_sat)
        implicit none
        real(real64) :: T_sat
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: p
        T_sat = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(p)) return
        if (p < self%region4%calc_psat(273.15d0) .or. p > 22.064d6) return
        T_sat = max(273.15d0, self%region4%calc_tsat(p))
        if (p == 22.064d6) T_sat = 647.096d0
    end function saturation_temperature_iapws97

    !> Batched saturation values and derivatives. density_dT is the total
    !> derivative along saturation; drho_dT and drho_dP are partials.
    module pure elemental subroutine calc_saturation_properties_iapws97(self, T_in, rho_vap, rho_liq, cp_vap, &
                                                                        cp_liq, drho_dT_vap, drho_dT_liq, drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, &
                                                                        dcp_dT_liq, latent_heat, dlatent_heat_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(out), optional :: rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, drho_dT_liq, &
                                               drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, latent_heat, &
                                               dlatent_heat_dT
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        ! Only heat-capacity derivatives require third-order EOS terms.
        call saturation_state(self, T_in, vap, liq, slope, &
                              present(dcp_dT_vap), present(dcp_dT_liq))
        if (present(rho_vap)) rho_vap = vap%rho
        if (present(rho_liq)) rho_liq = liq%rho
        if (present(cp_vap)) cp_vap = vap%cp
        if (present(cp_liq)) cp_liq = liq%cp
        if (present(drho_dT_vap)) drho_dT_vap = vap%rho_T
        if (present(drho_dT_liq)) drho_dT_liq = liq%rho_T
        if (present(drho_dP_vap)) drho_dP_vap = vap%rho_p
        if (present(drho_dP_liq)) drho_dP_liq = liq%rho_p
        if (present(density_dT_vap)) density_dT_vap = vap%rho_T + vap%rho_p * slope
        if (present(density_dT_liq)) density_dT_liq = liq%rho_T + liq%rho_p * slope
        if (present(dcp_dT_vap)) dcp_dT_vap = vap%cp_T + vap%cp_p * slope
        if (present(dcp_dT_liq)) dcp_dT_liq = liq%cp_T + liq%cp_p * slope
        if (present(latent_heat)) latent_heat = vap%h - liq%h
        if (present(dlatent_heat_dT)) dlatent_heat_dT = (vap%cp - liq%cp) + (vap%h_p - liq%h_p) * slope
    end subroutine

    module pure elemental subroutine calc_cp_derivatives_iapws97(self, T_in, p_in, dcp_dT, dcp_dp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT, dcp_dp
        call cp_derivatives(self, T_in, p_in, dcp_dT, dcp_dp)
    end subroutine

    !> Initialize the IAPWS-97 model instance.
    !> Initializes all sub-region handlers and auxiliary equations.
    module pure elemental subroutine initialize_type_iapws97(self)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(inout) :: self

        call self%auxiliary%initialize()
        call self%region1%initialize()
        call self%region2%initialize()
        call self%region3%initialize()
        call self%region4%initialize()
        call self%region5%initialize()

    end subroutine initialize_type_iapws97

    !> Determine the IAPWS-97 region for a given temperature and pressure.
    !> Checks validity limits and boundaries between regions (Saturation, B23).

    module pure elemental function get_region_iapws97(self, T_in, p_in) result(region_id)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Region ID result
        integer(int32) :: region_id

        real(real64) :: p_boundary

        ! 1. Global Validity Check
        ! Region 1, 2, 3: p <= 100 MPa
        ! Region 5: p <= 50 MPa -> Checked individually
        ! T >= 243.15 K
        region_id = IAPWS97_INVALID
        if (.not. ieee_is_finite(T_in) .or. .not. ieee_is_finite(p_in)) return
        if (p_in <= 0.0d0 .or. p_in > 100.0d6 .or. &
            T_in < 243.15d0 .or. T_in > 2273.15d0) then
            region_id = IAPWS97_INVALID
            return
        end if

        ! 2. Hierarchical determination by temperature (Based on Fig. 1)

        ! --- High Temperature (Region 5) ---
        ! 1073.15 K < T <= 2273.15 K
        if (T_in > 1073.15d0) then
            if (p_in <= 50.0d6) then
                region_id = IAPWS97_REGION_5
            else
                region_id = IAPWS97_INVALID ! Region 5 max P is 50MPa
            end if
            return
        end if

        ! --- Region 2 Only Area ---
        ! 863.15 K < T <= 1073.15 K
        ! This temp range is Region 2 for all pressures up to 100MPa
        if (T_in > 863.15d0) then
            region_id = IAPWS97_REGION_2
            return
        end if

        ! --- Boundary between Region 2 and Region 3 ---
        ! 623.15 K < T <= 863.15 K
        ! Boundary: B23 equation (Eq. 5)
        if (T_in > 623.15d0) then
            p_boundary = self%auxiliary%calc_p_boundary(T_in)

            if (p_in > p_boundary) then
                region_id = IAPWS97_REGION_3
            else
                region_id = IAPWS97_REGION_2
            end if
            return
        end if

        ! --- Boundary between Region 1 and Region 2 (Saturation Curve) ---
        ! 273.15 K <= T <= 623.15 K
        ! Boundary: Saturation pressure equation (Eq. 30)
        p_boundary = self%region4%calc_psat(T_in)

        if (p_in > p_boundary) then
            region_id = IAPWS97_REGION_1 ! Compressed liquid
        else
            region_id = IAPWS97_REGION_2 ! Superheated (or metastable) steam
        end if

    end function get_region_iapws97

    !> Calculate all thermodynamic properties for a given state.
    !> Dispatches calculation to the appropriate region based on \( T \) and \( P \).
    module pure elemental subroutine calc_properties_iapws97(self, T_in, p_in, property)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        type(type_iapws_property), intent(inout) :: property
        call self%evaluate_state(T_in, p_in, property)
    end subroutine calc_properties_iapws97

    module pure elemental subroutine evaluate_state_iapws97(self, T_in, p_in, property, guess_rho)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Output property structure
        type(type_iapws_property), intent(inout) :: property

        real(real64), intent(in), optional :: guess_rho
        real(real64) :: rho, nan

        property%region_id = self%get_region(T_in, p_in)
        select case (property%region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_2)
            call self%region2%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho, guess_rho)
            if (.not. ieee_is_finite(rho)) then
                nan = rho
                property = type_iapws_property(IAPWS97_REGION_3, nan, nan, nan, nan, nan, nan, nan, nan, &
                                               nan, nan, p_in, T_in, nan, nan, nan, nan)
                return
            end if
            call self%region3%calc_properties(T_in, rho, property)
        case (IAPWS97_REGION_5)
            call self%region5%calc_properties(T_in, p_in, property)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_properties_iapws97."
        end select

    end subroutine evaluate_state_iapws97

    !> Calculate specific volume \( \nu \).
    module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific volume [m^3/kg]
        real(real64), intent(inout) :: nu

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_nu(T_in, p_in, nu)
        case (IAPWS97_REGION_2)
            call self%region2%calc_nu(T_in, p_in, nu)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            nu = 1.0d0 / rho
        case (IAPWS97_REGION_5)
            call self%region5%calc_nu(T_in, p_in, nu)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_nu_iapws97."
        end select
    end subroutine calc_nu_iapws97

    !> Calculate density \( \rho \).
    module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Density [kg/m^3]
        real(real64), intent(inout) :: rho

        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_rho(T_in, p_in, rho)
        case (IAPWS97_REGION_2)
            call self%region2%calc_rho(T_in, p_in, rho)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
        case (IAPWS97_REGION_5)
            call self%region5%calc_rho(T_in, p_in, rho)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_rho_iapws97."
        end select

    end subroutine calc_rho_iapws97

    module pure elemental subroutine calc_drho_dT_iapws97(self, T_in, p_in, drho_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dT

        real(real64) :: rho
        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_drho_dT(T_in, rho, drho_dT)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dT(T_in, p_in, drho_dT)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_drho_dT_iapws97."
        end select

    end subroutine calc_drho_dT_iapws97

    module pure elemental subroutine calc_drho_dp_iapws97(self, T_in, p_in, drho_dp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dp

        real(real64) :: rho
        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_drho_dp(T_in, rho, drho_dp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dp(T_in, p_in, drho_dp)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_drho_dp_iapws97."
        end select

    end subroutine calc_drho_dp_iapws97

    !> Calculate specific internal energy \( u \).
    module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific internal energy [J/kg]
        real(real64), intent(inout) :: u

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_u(T_in, p_in, u)
        case (IAPWS97_REGION_2)
            call self%region2%calc_u(T_in, p_in, u)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_u(T_in, rho, u)
        case (IAPWS97_REGION_5)
            call self%region5%calc_u(T_in, p_in, u)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_u_iapws97."
        end select

    end subroutine calc_u_iapws97

    !> Calculate specific enthalpy \( h \).
    module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific enthalpy [J/kg]
        real(real64), intent(inout) :: h

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_h(T_in, p_in, h)
        case (IAPWS97_REGION_2)
            call self%region2%calc_h(T_in, p_in, h)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_h(T_in, rho, h)
        case (IAPWS97_REGION_5)
            call self%region5%calc_h(T_in, p_in, h)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_h_iapws97."
        end select

    end subroutine calc_h_iapws97

    !> Calculate specific entropy \( s \).
    module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific entropy [J/(kg K)]
        real(real64), intent(inout) :: s

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_s(T_in, p_in, s)
        case (IAPWS97_REGION_2)
            call self%region2%calc_s(T_in, p_in, s)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_s(T_in, rho, s)
        case (IAPWS97_REGION_5)
            call self%region5%calc_s(T_in, p_in, s)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_s_iapws97."
        end select

    end subroutine calc_s_iapws97

    !> Calculate specific isobaric heat capacity \( c_p \).
    module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific isobaric heat capacity [J/(kg K)]
        real(real64), intent(inout) :: cp

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)
        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_cp(T_in, p_in, cp)
        case (IAPWS97_REGION_2)
            call self%region2%calc_cp(T_in, p_in, cp)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_cp(T_in, rho, cp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cp(T_in, p_in, cp)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_cp_iapws97."
        end select

    end subroutine calc_cp_iapws97

    !> Calculate specific isochoric heat capacity \( c_v \).
    module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific isochoric heat capacity [J/(kg K)]
        real(real64), intent(inout) :: cv

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)
        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_cv(T_in, p_in, cv)
        case (IAPWS97_REGION_2)
            call self%region2%calc_cv(T_in, p_in, cv)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_cv(T_in, rho, cv)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cv(T_in, p_in, cv)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_cv_iapws97."
        end select

    end subroutine calc_cv_iapws97

    !> Calculate speed of sound \( w \).
    module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Speed of sound [m/s]
        real(real64), intent(inout) :: w

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_w(T_in, p_in, w)
        case (IAPWS97_REGION_2)
            call self%region2%calc_w(T_in, p_in, w)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_w(T_in, rho, w)
        case (IAPWS97_REGION_5)
            call self%region5%calc_w(T_in, p_in, w)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_w_iapws97."
        end select

    end subroutine calc_w_iapws97

    !> cp temperature derivative at fixed pressure [J/(kg K^2)].
    module pure elemental subroutine calc_dcp_dT_iapws97(self, T_in, p_in, dcp_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT
        real(real64) :: unused
        call cp_derivatives(self, T_in, p_in, dcp_dT, unused)
    end subroutine

    !> cp pressure derivative at fixed temperature [J/(kg K Pa)].
    module pure elemental subroutine calc_dcp_dp_iapws97(self, T_in, p_in, dcp_dp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dp
        real(real64) :: unused
        call cp_derivatives(self, T_in, p_in, unused, dcp_dp)
    end subroutine

    pure elemental subroutine cp_derivatives(self, T, p, ct, cp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T, p
        real(real64), intent(out) :: ct, cp
        real(real64) :: rho
        ct = ieee_value(0d0, ieee_quiet_nan); cp = ct
        if (.not. ieee_is_finite(T) .or. .not. ieee_is_finite(p)) return
        if (T == 647.096d0 .and. abs(p - 22.064d6) < 1d-3) return
        select case (self%get_region(T, p))
        case (1)
            call self%region1%calc_cp_derivatives(T, p, ct, cp)
        case (2)
            call self%region2%calc_cp_derivatives(T, p, ct, cp)
        case (3)
            call region3_density(self, T, p, rho)
            call self%region3%calc_cp_derivatives(T, rho, ct, cp)
        case (5)
            call self%region5%calc_cp_derivatives(T, p, ct, cp)
        end select
    end subroutine

    !> Safeguarded region-3 solve shared by values and ordinary derivatives.
    pure subroutine evaluate_region3_pressure(delta, coef, reduced_p, reduced_dp)
        implicit none
        real(real64), intent(in) :: delta, coef(0:11)
        real(real64), intent(out) :: reduced_p, reduced_dp
        integer :: i

        reduced_p = coef(11)
        reduced_dp = 12d0 * coef(11)
        do i = 10, 0, -1
            reduced_p = coef(i) + delta * reduced_p
            reduced_dp = (i + 1) * coef(i) + delta * reduced_dp
        end do
        reduced_p = delta * reduced_p
    end subroutine evaluate_region3_pressure

    pure elemental subroutine region3_density(self, T, p, rho, guess_rho)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T, p
        real(real64), intent(out) :: rho
        real(real64), intent(in), optional :: guess_rho
        real(real64) :: lo, hi, mid, pm, slope, trial, estimate, Ts
        real(real64) :: coefficients(0:11), rt, reduced_p, reduced_dp, delta, ddp, correction
        integer :: k
        logical :: vapor
        integer :: i
        vapor = .false.
        if (T < self%region3%T_c) vapor = p <= self%region4%calc_psat(T)
        ! SR5-05 seeds bounded Newton refinement against the forward EOS.
        ! Keep the bracketed solver at the critical point and as a fallback.
        if (abs(T - self%region3%T_c) > 0.1d0 .or. abs(p - 22.064d6) > 1d5) then
            ! For region selection the seed only compares T with Tsat.
            ! Reuse the side already selected by the saturation pressure.
            Ts = merge(T - 1d0, T + 1d0, vapor .or. T >= self%region3%T_c)
            estimate = 0d0
            if (present(guess_rho)) then
                if (ieee_is_finite(guess_rho)) then
                    if (guess_rho > 1d0 .and. guess_rho < 800d0) then
                        if (T >= self%region3%T_c .or. ((guess_rho < self%region3%rho_c) .eqv. vapor)) &
                            estimate = 1d0 / guess_rho
                    end if
                end if
            end if
            if (estimate <= 0d0) estimate = volume_region3(T, p / 1d6, Ts)
            if (ieee_is_finite(estimate)) then
                if (estimate > 0d0) then
                    mid = 1d0 / estimate
                    rt = self%region3%R * T
                    call prepare_region3_pressure(self%region3%T_c / T, coefficients)
                    do i = 1, 8
                        if (mid < 1d0 .or. mid > 800d0) exit
                        if (T < self%region3%T_c) then
                            if ((mid < self%region3%rho_c) .neqv. vapor) exit
                        end if
                        delta = mid / self%region3%rho_c
                        call evaluate_region3_pressure(delta, coefficients, reduced_p, reduced_dp)
                        pm = reduced_p * (self%region3%rho_c * rt)
                        slope = reduced_dp * rt
                        if (slope <= 0d0) exit
                        if (abs(pm - p) <= max(2d-7, 16d0 * epsilon(mid) * mid * slope)) then
                            ! Roundoff in the grouped coefficients is magnified
                            ! in the dense liquid and on small-slope isotherms.
                            if (mid > 400d0 .and. (T < 625d0 .or. slope < 1.1d4)) then
                                call self%region3%calc_pressure_pair(T, mid, pm, slope)
                                if (slope <= 0d0 .or. abs(pm - p) > 1d-3) exit
                                mid = mid - (pm - p) / slope
                            end if
                            rho = mid
                            return
                        end if
                        correction = (pm - p) / slope
                        if (i == 1) then
                            ! Use a guarded Halley correction for the initial fit.
                            ddp = 132d0 * coefficients(11)
                            do k = 10, 1, -1
                                ddp = k * (k + 1) * coefficients(k) + delta * ddp
                            end do
                            trial = 0.5d0 * correction * ddp / (self%region3%rho_c * reduced_dp)
                            if (abs(trial) < 0.125d0) correction = correction / (1d0 - trial)
                        end if
                        if (abs(correction) > 0.05d0 * mid) exit
                        mid = mid - correction
                    end do
                end if
            end if
        end if
        if (T < self%region3%T_c) then
            call saturation_root(self, T, p, vapor, rho)
            return
        end if
        lo = 1d0; hi = 800d0; mid = (lo + hi) / 2d0
        do i = 1, 100
            call self%region3%calc_pressure_pair(T, mid, pm, slope)
            if (pm > p) then
                hi = mid
            else
                lo = mid
            end if
            if (hi - lo < 2d-11) exit
            trial = (lo + hi) / 2d0
            if (slope > 0d0) then
                if (mid - (pm - p) / slope > lo .and. mid - (pm - p) / slope < hi) trial = mid - (pm - p) / slope
            end if
            mid = trial
        end do
        rho = (lo + hi) / 2d0
        call self%region3%calc_pressure_pair(T, rho, pm, slope)
        if (abs(pm - p) > 0.01d0) rho = ieee_value(0d0, ieee_quiet_nan)
    end subroutine

    !> Stable phase root, bracketed by spinodals to prevent phase switching.
    pure elemental subroutine saturation_root(self, T, p, vapor, rho)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T, p
        logical, intent(in) :: vapor
        real(real64), intent(out) :: rho
        real(real64) :: lo, hi, mid, slope, plo, phi, pm, trial
        integer :: i
        rho = ieee_value(0d0, ieee_quiet_nan)
        ! Locate the spinodal on the requested side of the critical density.
        if (vapor) then
            lo = 1d0; hi = self%region3%rho_c
        else
            lo = self%region3%rho_c; hi = 800d0
        end if
        do i = 1, 70
            mid = (lo + hi) / 2d0
            call self%region3%calc_pressure_pair(T, mid, pm, slope)
            if ((slope > 0d0) .eqv. vapor) then
                lo = mid
            else
                hi = mid
            end if
        end do
        if (vapor) then
            hi = (lo + hi) / 2d0; lo = 1d0
        else
            lo = (lo + hi) / 2d0; hi = 800d0
        end if
        call self%region3%calc_pressure_pair(T, lo, plo, slope)
        call self%region3%calc_pressure_pair(T, hi, phi, slope)
        if (p < plo .or. p > phi) return
        mid = (lo + hi) / 2d0
        do i = 1, 100
            call self%region3%calc_pressure_pair(T, mid, pm, slope)
            if (pm > p) then
                hi = mid
            else
                lo = mid
            end if
            if (hi - lo < 2d-11) exit
            trial = (lo + hi) / 2d0
            if (slope > 0d0) then
                if (mid - (pm - p) / slope > lo .and. mid - (pm - p) / slope < hi) trial = mid - (pm - p) / slope
            end if
            mid = trial
        end do
        rho = (lo + hi) / 2d0
        call self%region3%calc_pressure_pair(T, rho, pm, slope)
        if (abs(pm - p) > 0.01d0 .or. slope <= 0d0) rho = ieee_value(0d0, ieee_quiet_nan)
    end subroutine

    pure elemental subroutine saturation_state(self, T, vap, liq, slope, need_cp_derivatives_vap, need_cp_derivatives_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T
        type(saturation_phase), intent(out) :: vap, liq
        real(real64), intent(out) :: slope
        logical, intent(in) :: need_cp_derivatives_vap, need_cp_derivatives_liq
        type(type_iapws_gibbs_coefficient) :: g
        real(real64) :: p, nan
        nan = ieee_value(0d0, ieee_quiet_nan)
        vap = saturation_phase(nan, nan, nan, nan, nan, nan, nan, nan); liq = vap; slope = nan
        if (.not. ieee_is_finite(T)) return
        if (T < 273.15d0 .or. T > 647.096d0) return
        if (T == 647.096d0) then
            vap%rho = 322d0; liq%rho = 322d0
            call self%region3%calc_h(T, 322d0, vap%h)
            liq%h = vap%h
            return
        end if
        p = self%region4%calc_psat(T)
        slope = self%region4%calc_dpsat_dT(T)
        if (T <= 623.15d0) then
            if (need_cp_derivatives_vap) then
                call self%region2%calc_gamma_cp(T, p, g, vap%cp_T, vap%cp_p)
            else
                call self%region2%calc_gamma(T, p, g)
            end if
            call saturation_gibbs(T, g, vap)
            if (need_cp_derivatives_liq) then
                call self%region1%calc_gamma_cp(T, p, g, liq%cp_T, liq%cp_p)
            else
                call self%region1%calc_gamma(T, p, g)
            end if
            call saturation_gibbs(T, g, liq)
        else
            call saturation_root(self, T, p, .true., vap%rho)
            call saturation_root(self, T, p, .false., liq%rho)
            if (.not. ieee_is_finite(vap%rho) .or. .not. ieee_is_finite(liq%rho)) then
                ! The independent IF97 fits can lose a stable root extremely
                ! close to Tc. Reject the whole pair instead of leaking one phase.
                vap%rho = nan
                liq%rho = nan
                return
            end if
            call saturation_helmholtz(self%region3, T, vap, need_cp_derivatives_vap)
            call saturation_helmholtz(self%region3, T, liq, need_cp_derivatives_liq)
        end if
    end subroutine

    pure elemental subroutine saturation_gibbs(T, g, phase)
        implicit none
        real(real64), intent(in) :: T
        type(saturation_phase), intent(inout) :: phase
        type(type_iapws_gibbs_coefficient), intent(in) :: g
        phase%rho = 1d0 / g%g_p
        phase%cp = -T * g%g_tt
        phase%h = g%g - T * g%g_t
        phase%rho_T = -g%g_pt / g%g_p**2
        phase%rho_p = -g%g_pp / g%g_p**2
        phase%h_p = g%g_p - T * g%g_pt
    end subroutine

    pure elemental subroutine saturation_helmholtz(region, T, phase, need_cp_derivatives)
        implicit none
        class(type_iapws97_region3), intent(in) :: region
        real(real64), intent(in) :: T
        type(saturation_phase), intent(inout) :: phase
        logical, intent(in) :: need_cp_derivatives
        type(type_iapws_helmholtz_property) :: f
        real(real64) :: pt, pr
        call region%calc_phi(region%T_c / T, phase%rho / region%rho_c, f)
        call region%calc_cp(T, phase%rho, phase%cp, f)
        call region%calc_h(T, phase%rho, phase%h, f)
        call region%calc_p_T(T, phase%rho, pt, f)
        call region%calc_p_rho(T, phase%rho, pr, f)
        phase%rho_T = -pt / pr
        phase%rho_p = 1d0 / pr
        ! h_p|T = v - T*v_T|p, valid also for the Helmholtz EOS.
        phase%h_p = 1d0 / phase%rho + T * phase%rho_T / phase%rho**2
        if (need_cp_derivatives) call region%calc_cp_derivatives(T, phase%rho, phase%cp_T, phase%cp_p, f)
    end subroutine

    module pure elemental subroutine calc_latent_heat_iapws97(self, latent_heat, T_in, p_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(inout) :: latent_heat
        real(real64), intent(in), optional :: T_in, p_in
        type(saturation_phase) :: vap, liq
        real(real64) :: T, slope
        latent_heat = ieee_value(0d0, ieee_quiet_nan)
        if (present(T_in) .eqv. present(p_in)) return
        if (present(T_in)) then
            T = T_in
        else
            if (.not. ieee_is_finite(p_in)) return
            if (p_in < self%region4%calc_psat(273.15d0) .or. p_in > 22.064d6) return
            T = self%region4%calc_tsat(p_in)
            if (p_in == 22.064d6) T = 647.096d0
        end if
        call saturation_state(self, T, vap, liq, slope, .false., .false.)
        latent_heat = vap%h - liq%h
    end subroutine

    module pure elemental subroutine calc_latent_heat_dT_iapws97(self, T_in, dlatent_heat_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(out) :: dlatent_heat_dT
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        dlatent_heat_dT = (vap%cp - liq%cp) + (vap%h_p - liq%h_p) * slope
    end subroutine

    module pure elemental subroutine calc_saturation_density_iapws97(self, T_in, rho_vap, rho_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: rho_vap, rho_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        if (present(rho_vap)) rho_vap = vap%rho
        if (present(rho_liq)) rho_liq = liq%rho
    end subroutine

    module pure elemental subroutine calc_saturation_cp_iapws97(self, T_in, cp_vap, cp_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: cp_vap, cp_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        if (present(cp_vap)) cp_vap = vap%cp
        if (present(cp_liq)) cp_liq = liq%cp
    end subroutine

    module pure elemental subroutine calc_saturation_drho_dT_iapws97(self, T_in, drho_dT_vap, drho_dT_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: drho_dT_vap, drho_dT_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        if (present(drho_dT_vap)) drho_dT_vap = vap%rho_T
        if (present(drho_dT_liq)) drho_dT_liq = liq%rho_T
    end subroutine

    module pure elemental subroutine calc_saturation_drho_dP_iapws97(self, T_in, drho_dP_vap, drho_dP_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: drho_dP_vap, drho_dP_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        if (present(drho_dP_vap)) drho_dP_vap = vap%rho_p
        if (present(drho_dP_liq)) drho_dP_liq = liq%rho_p
    end subroutine

    module pure elemental subroutine calc_saturation_dcp_dT_iapws97(self, T_in, dcp_dT_vap, dcp_dT_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(out), optional :: dcp_dT_vap, dcp_dT_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, present(dcp_dT_vap), present(dcp_dT_liq))
        if (present(dcp_dT_vap)) dcp_dT_vap = vap%cp_T + vap%cp_p * slope
        if (present(dcp_dT_liq)) dcp_dT_liq = liq%cp_T + liq%cp_p * slope
    end subroutine

    module pure elemental subroutine calc_saturation_density_dT_iapws97(self, T_in, drho_dT_vap, drho_dT_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(out), optional :: drho_dT_vap, drho_dT_liq
        type(saturation_phase) :: vap, liq
        real(real64) :: slope
        call saturation_state(self, T_in, vap, liq, slope, .false., .false.)
        if (present(drho_dT_vap)) drho_dT_vap = vap%rho_T + vap%rho_p * slope
        if (present(drho_dT_liq)) drho_dT_liq = liq%rho_T + liq%rho_p * slope
    end subroutine

end submodule iapws97_base
