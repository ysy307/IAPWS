submodule(module_iapws97) iapws97_base
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
contains

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
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Output property structure
        type(type_iapws_property), intent(inout) :: property

        property%region_id = self%get_region(T_in, p_in)
        select case (property%region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_2)
            call self%region2%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_3)
            call self%region3%calc_rho(T_in, p_in, property%rho)
            call self%region3%calc_properties(T_in, property%rho, property)
        case (IAPWS97_REGION_5)
            call self%region5%calc_properties(T_in, p_in, property)
        end select

    end subroutine calc_properties_iapws97

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
            call self%region3%calc_rho(T_in, p_in, rho)
            nu = 1.0d0 / rho
        case (IAPWS97_REGION_5)
            call self%region5%calc_nu(T_in, p_in, nu)
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
            call self%region3%calc_rho(T_in, p_in, rho)
        case (IAPWS97_REGION_5)
            call self%region5%calc_rho(T_in, p_in, rho)
        end select

    end subroutine calc_rho_iapws97

    module pure elemental subroutine calc_drho_dT_iapws97(self, T_in, p_in, drho_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dT

        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (self%get_region(T_in, p_in))
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_3)
            call self%region3%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dT(T_in, p_in, drho_dT)
        end select

    end subroutine calc_drho_dT_iapws97

    module pure elemental subroutine calc_drho_dp_iapws97(self, T_in, p_in, drho_dp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dp

        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (self%get_region(T_in, p_in))
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_3)
            call self%region3%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dp(T_in, p_in, drho_dp)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_u(T_in, rho, u)
        case (IAPWS97_REGION_5)
            call self%region5%calc_u(T_in, p_in, u)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_h(T_in, rho, h)
        case (IAPWS97_REGION_5)
            call self%region5%calc_h(T_in, p_in, h)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_s(T_in, rho, s)
        case (IAPWS97_REGION_5)
            call self%region5%calc_s(T_in, p_in, s)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_cp(T_in, rho, cp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cp(T_in, p_in, cp)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_cv(T_in, rho, cv)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cv(T_in, p_in, cv)
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
            call self%region3%calc_rho(T_in, p_in, rho)
            call self%region3%calc_w(T_in, rho, w)
        case (IAPWS97_REGION_5)
            call self%region5%calc_w(T_in, p_in, w)
        end select

    end subroutine calc_w_iapws97

    !> 蒸発潜熱の計算 (Implementation)
    module pure elemental subroutine calc_latent_heat_iapws97(self, latent_heat, T_in, p_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(inout) :: latent_heat
        real(real64), intent(in), optional :: T_in
        real(real64), intent(in), optional :: p_in

        real(real64) :: T_sat, P_sat
        real(real64) :: h_liq, h_vap

        ! 1. 入力に基づき飽和点 (T, P) を特定
        if (present(T_in) .and. .not. present(p_in)) then
            ! 温度指定 -> 飽和圧力を計算 (Region 4)
            T_sat = T_in
            P_sat = self%region4%calc_psat(T_sat)

        else if (.not. present(T_in) .and. present(p_in)) then
            ! 圧力指定 -> 飽和温度を計算 (Region 4)
            P_sat = p_in
            T_sat = self%region4%calc_tsat(P_sat)

        else
            ! エラー: 両方指定、または両方なし
            ! 必要に応じて error stop 等を入れてください
            latent_heat = 0.0_real64
            return
        end if

        ! 2. 飽和液体エンタルピー h' (Region 1)
        ! Region 1 は abst_iapws_gibbs を継承しているため calc_h が使えます
        call self%region1%calc_h(T_sat, P_sat, h_liq)

        ! 3. 飽和蒸気エンタルピー h'' (Region 2)
        ! Region 2 も abst_iapws_gibbs を継承しているため calc_h が使えます
        call self%region2%calc_h(T_sat, P_sat, h_vap)

        ! 4. 潜熱 L = h'' - h'
        latent_heat = h_vap - h_liq

    end subroutine calc_latent_heat_iapws97

    !> 飽和蒸気密度と飽和液体密度を計算 (T入力)
    module pure elemental subroutine calc_saturation_density_iapws97(self, T_in, rho_vap, rho_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: rho_vap ! 飽和蒸気密度 rho''
        real(real64), intent(inout), optional :: rho_liq ! 飽和液体密度 rho'

        real(real64) :: P_sat

        ! 1. 飽和圧力 P_sat を計算 (Region 4)
        P_sat = self%region4%calc_psat(T_in)

        ! 2. 飽和蒸気密度 (Region 2)
        if (present(rho_vap)) then
            ! Region 2 の calc_rho を呼ぶ
            call self%region2%calc_rho(T_in, P_sat, rho_vap)
        end if

        ! 3. 飽和液体密度 (Region 1)
        if (present(rho_liq)) then
            ! Region 1 の calc_rho を呼ぶ
            call self%region1%calc_rho(T_in, P_sat, rho_liq)
        end if

    end subroutine calc_saturation_density_iapws97

    module pure elemental subroutine calc_saturation_cp(self, T_in, cp_vap, cp_liq)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(inout), optional :: cp_vap ! 飽和蒸気比熱 cp''
        real(real64), intent(inout), optional :: cp_liq ! 飽和液体比熱 cp'

        real(real64) :: P_sat

        ! 1. 飽和圧力 P_sat を計算 (Region 4)
        P_sat = self%region4%calc_psat(T_in)

        ! 2. 飽和蒸気比熱 (Region 2)
        if (present(cp_vap)) then
            ! Region 2 の calc_cp を呼ぶ
            call self%region2%calc_cp(T_in, P_sat, cp_vap)
        end if

        ! 3. 飽和液体比熱 (Region 1)
        if (present(cp_liq)) then
            ! Region 1 の calc_cp を呼ぶ
            call self%region1%calc_cp(T_in, P_sat, cp_liq)
        end if

    end subroutine calc_saturation_cp

end submodule iapws97_base
