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
            call self%region3%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_5)
            call self%region5%calc_properties(T_in, p_in, property)
        end select

    end subroutine calc_properties_iapws97

    module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: nu
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_nu_iapws97

    module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: rho
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_rho_iapws97

    module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: u
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_u_iapws97

    module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: h
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_h_iapws97

    module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: s
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_s_iapws97

    module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: cp
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_cp_iapws97

    module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: cv
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_cv_iapws97

    module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w, prop_in)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: w
        type(type_iapws_property), intent(inout), optional :: prop_in

    end subroutine calc_w_iapws97

end submodule iapws97_base
