submodule(module_iapws97) iapws97_base
    implicit none
    integer(int32), parameter :: IAPWS97_INVALID = -1
    integer(int32), parameter :: IAPWS97_REGION_1 = 1
    integer(int32), parameter :: IAPWS97_REGION_2 = 2
    integer(int32), parameter :: IAPWS97_REGION_3 = 3
    integer(int32), parameter :: IAPWS97_REGION_5 = 5

    real(real64), parameter :: IAPWS97_LIMIT_T_MAX = 2273.15d0
    real(real64), parameter :: IAPWS97_LIMIT_P_MAX = 100.0d6

    real(real64), parameter :: IAPWS97_R5_T_MIN = 1073.15d0
    real(real64), parameter :: IAPWS97_R5_P_MAX = 50.0d6
contains

    module pure elemental subroutine initialize_type_iapws97(self)
        implicit none
        class(type_iapws97), intent(inout) :: self

        call self%auxiliary%initialize()
        call self%region1%initialize()
        call self%region2%initialize()
        call self%region3%initialize()
        call self%region4%initialize()
        call self%region5%initialize()

    end subroutine initialize_type_iapws97

    module pure elemental function get_region_iapws97(self, T_in, p_in) result(region_id)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        integer(int32) :: region_id

        real(real64) :: p_boundary

        ! 1. 全体の有効範囲チェック (Global Validity)
        ! Region 1, 2, 3: p <= 100 MPa [cite: 119]
        ! Region 5: p <= 50 MPa [cite: 1028] -> 個別に判定
        ! T >= 243.15 K [cite: 117]
        if (p_in <= 0.0d0 .or. p_in > 100.0d6 .or. &
            T_in < 243.15d0 .or. T_in > 2273.15d0) then
            region_id = IAPWS97_INVALID
            return
        end if

        ! 2. 温度による階層的な判定 (Fig. 1に基づく分岐)

        ! --- 高温領域 (Region 5) ---
        ! 1073.15 K < T <= 2273.15 K [cite: 1027]
        if (T_in > 1073.15d0) then
            if (p_in <= 50.0d6) then
                region_id = IAPWS97_REGION_5
            else
                region_id = IAPWS97_INVALID ! Region 5の上限は50MPa
            end if
            return
        end if

        ! --- Region 2 の単独領域 ---
        ! 863.15 K < T <= 1073.15 K [cite: 482]
        ! この温度帯は全圧域(100MPaまで)でRegion 2確定
        if (T_in > 863.15d0) then
            region_id = IAPWS97_REGION_2
            return
        end if

        ! --- Region 2 と Region 3 の境界 ---
        ! 623.15 K < T <= 863.15 K [cite: 481, 901]
        ! 境界線: B23式 (Eq. 5) [cite: 172]
        if (T_in > 623.15d0) then
            ! B23式 (2次多項式) は計算が軽いのでここで呼ぶ
            p_boundary = self%auxiliary%calc_p_boundary(T_in)

            if (p_in > p_boundary) then
                region_id = IAPWS97_REGION_3
            else
                region_id = IAPWS97_REGION_2
            end if
            return
        end if

        ! --- Region 1 と Region 2 の境界 (飽和曲線) ---
        ! 273.15 K <= T <= 623.15 K [cite: 353, 480]
        ! 境界線: 飽和蒸気圧式 (Eq. 30) [cite: 921]
        ! ※この計算が最も重いため、最後に判定する

        p_boundary = self%region4%calc_psat(T_in)

        if (p_in > p_boundary) then
            region_id = IAPWS97_REGION_1 ! 圧縮水
        else
            region_id = IAPWS97_REGION_2 ! 過熱蒸気 (または準安定蒸気)
        end if

    end function get_region_iapws97

    module pure elemental subroutine calc_properties_iapws97(self, T_in, p_in, property)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

    module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
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

end submodule iapws97_base
