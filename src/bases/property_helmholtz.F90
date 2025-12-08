submodule(module_iapws) iapws_property_helmholtz
    implicit none
contains

    module pure elemental subroutine calc_properties_helmholtz(self, T_in, rho_in, property)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        type(type_iapws_property), intent(inout) :: property

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        property%T = T_in
        property%rho = rho_in

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        call self%calc_phi(tau, delta, props)

        call self%calc_p(T_in, rho_in, property%p, props)
        call self%calc_u(T_in, rho_in, property%u, props)
        call self%calc_s(T_in, rho_in, property%s, props)
        call self%calc_h(T_in, rho_in, property%h, props)
        call self%calc_cp(T_in, rho_in, property%cp, props)
        call self%calc_cv(T_in, rho_in, property%cv, props)
        call self%calc_w(T_in, rho_in, property%w, props)

    end subroutine calc_properties_helmholtz

    module pure elemental subroutine calc_rho_helmholtz(self, T_in, P_in, rho)
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in, P_in
        real(real64), intent(inout) :: rho

        integer(int32) :: iter
        real(real64) :: p, dp, delta, tau
        real(real64) :: drho, rho_old
        real(real64) :: p_diff, threshold
        type(type_iapws_helmholtz_property) :: props

        ! --- 1. 初期値推測 (変更なし) ---
        if (P_in > 22.1d6) then
            rho = 500.0d0
        elseif (P_in > 10.0d6) then
            rho = 322.0d0
        else
            rho = P_in / (self%R * T_in)
        end if

        tau = self%T_c / T_in

        ! --- 2. 限界まで追い込む反復計算 ---
        do iter = 1, 100
            delta = rho / self%rho_c
            call self%calc_phi(tau, delta, props)
            call self%calc_p(T_in, rho, p, props)
            call self%calc_p_rho(T_in, rho, dp, props)

            ! ガード: 極端に小さな傾きによる発散防止
            if (dp < 1.0d-4) dp = 1.0d-4

            ! ニュートンステップ計算
            p_diff = p - P_in
            drho = -p_diff / dp
            rho_old = rho

            ! 現在の誤差の程度を確認
            ! (目標圧力に対する相対誤差)
            threshold = abs(p_diff) / max(abs(P_in), tiny(1.0d0))

            ! --- [重要] 2段階スイッチング ---
            if (threshold > 1.0d-2) then
                ! [A] まだ遠い場合 (誤差 > 1%):
                ! 暴走を防ぐため、変化量を制限する (ダンピング)
                if (drho > 0.2d0 * rho) drho = 0.2d0 * rho
                if (drho < -0.2d0 * rho) drho = -0.2d0 * rho
            end if

            rho = rho + drho

            ! 負の密度ガード
            if (rho <= 0.0d0) rho = 0.5d0 * rho_old

            ! --- [究極の終了判定] ---
            ! 1. 密度の変化量がマシンイプシロン(約2e-16)レベルになったか？
            ! 2. または、圧力の相対誤差がマシンイプシロンレベルになったか？
            ! どちらかを満たせば、「これ以上計算しても変わらない」ので終了。
            if (abs(drho) <= epsilon(rho) * rho) exit
            if (abs(p - P_in) <= epsilon(P_in) * P_in) exit

        end do
    end subroutine calc_rho_helmholtz

    module pure elemental subroutine calc_p_helmholtz(self, T_in, rho_in, p, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        p = rho_in * self%R * T_in * delta * props%f_d

    end subroutine calc_p_helmholtz

    module pure elemental subroutine calc_p_rho_helmholtz(self, T_in, rho_in, p_rho, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p_rho
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props
        real(real64) :: R, RT

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c
        R = self%R
        RT = R * T_in

        ! φ, φ_d, φ_dd の取得
        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        ! p_rho = R * T_in * (delta * props%f_d + delta * props%f_d + delta**2 * props%f_dd) / self%rho_c
        p_rho = R * T_in * (2.0d0 * delta * props%f_d + delta**2 * props%f_dd)

    end subroutine calc_p_rho_helmholtz

    module pure elemental subroutine calc_p_T_helmholtz(self, T_in, rho_in, p_T, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p_T
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props
        real(real64) :: R

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c
        R = self%R

        ! φ, φ_d, φ_dd, φ_dt の取得
        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        ! ∂p/∂T  (rho constant)
        ! p_T = rho R delta ( phi_d - tau * phi_dt )
        p_T = rho_in * R * delta * (props%f_d - tau * props%f_dt)

    end subroutine calc_p_T_helmholtz

    module pure elemental subroutine calc_u_helmholtz(self, T_in, rho_in, u, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: u
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        u = self%R * T_in * tau * props%f_t

    end subroutine calc_u_helmholtz

    module pure elemental subroutine calc_s_helmholtz(self, T_in, rho_in, s, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: s
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        s = self%R * (tau * props%f_t - props%f)

    end subroutine calc_s_helmholtz

    module pure elemental subroutine calc_h_helmholtz(self, T_in, rho_in, h, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: h
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        h = self%R * T_in * (tau * props%f_t + delta * props%f_d)

    end subroutine calc_h_helmholtz

    module pure elemental subroutine calc_cv_helmholtz(self, T_in, rho_in, cv, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: cv
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        cv = self%R * (-tau**2 * props%f_tt)

    end subroutine calc_cv_helmholtz

    module pure elemental subroutine calc_cp_helmholtz(self, T_in, rho_in, cp, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: cp
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        cp = self%R * (-tau**2 * props%f_tt + &
                       (delta * props%f_d - delta * tau * props%f_dt)**2 / &
                       (2.0d0 * delta * props%f_d + delta**2 * props%f_dd))

    end subroutine calc_cp_helmholtz

    module pure elemental subroutine calc_w_helmholtz(self, T_in, rho_in, w, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: w
        type(type_iapws_helmholtz_property), intent(in), optional :: prop_in

        real(real64) :: tau, delta
        real(real64) :: w_sq
        type(type_iapws_helmholtz_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        w_sq = self%R * T_in * &
               (2.0d0 * delta * props%f_d + delta**2 * props%f_dd - &
                ((delta * props%f_d - delta * tau * props%f_dt)**2) / &
                (tau**2 * props%f_tt))

        w = sqrt(max(w_sq, 0.0d0))
    end subroutine calc_w_helmholtz

end submodule iapws_property_helmholtz
