submodule(iapws) iapws_property_helmholtz
    implicit none
contains

    module pure elemental subroutine calc_properties_helmholtz(self, T_in, rho_in, property)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        type(type_iapws_property), intent(inout) :: property

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

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
        type(type_iapws_phi_property) :: props

        rho = p_in / (self%R * T_in)

        tau = self%T_c / T_in

        do iter = 1, 100
            delta = rho / self%rho_c
            call self%calc_phi(tau, delta, props)
            call self%calc_p(T_in, rho, p, props)
            call self%calc_p_rho(T_in, rho, dp, props)

            drho = -(p - p_in) / dp
            rho_old = rho
            rho = rho + drho

            if (abs(drho) < 1.0d-12 * rho_old) exit
        end do
    end subroutine calc_rho_helmholtz

    module pure elemental subroutine calc_p_helmholtz(self, T_in, rho_in, p, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        p = rho_in * self%R * T_in * delta * props%phi_d

    end subroutine calc_p_helmholtz

    module pure elemental subroutine calc_p_rho_helmholtz(self, T_in, rho_in, p_rho, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p_rho
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props
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

        ! dp/drho の解析式
        p_rho = RT * (delta * props%phi_d &
                      + rho_in * ((1.0d0 / self%rho_c) * props%phi_d &
                                  + delta * props%phi_dd / self%rho_c))

    end subroutine calc_p_rho_helmholtz

    module pure elemental subroutine calc_p_T_helmholtz(self, T_in, rho_in, p_T, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: p_T
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props
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
        p_T = rho_in * R * delta * (props%phi_d - tau * props%phi_dt)

    end subroutine calc_p_T_helmholtz

    module pure elemental subroutine calc_u_helmholtz(self, T_in, rho_in, u, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: u
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        u = self%R * T_in * tau * props%phi_t

    end subroutine calc_u_helmholtz

    module pure elemental subroutine calc_s_helmholtz(self, T_in, rho_in, s, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: s
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        s = self%R * (tau * props%phi_t - props%phi)

    end subroutine calc_s_helmholtz

    module pure elemental subroutine calc_h_helmholtz(self, T_in, rho_in, h, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: h
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        h = self%R * T_in * (tau * props%phi_t + delta * props%phi_d)

    end subroutine calc_h_helmholtz

    module pure elemental subroutine calc_cv_helmholtz(self, T_in, rho_in, cv, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: cv
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        cv = self%R * (-tau**2 * props%phi_tt)

    end subroutine calc_cv_helmholtz

    module pure elemental subroutine calc_cp_helmholtz(self, T_in, rho_in, cp, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: cp
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        cp = self%R * (-tau**2 * props%phi_tt + &
                       (delta * props%phi_d - delta * tau * props%phi_dt)**2 / &
                       (2.0d0 * delta * props%phi_d + delta**2 * props%phi_dd))

    end subroutine calc_cp_helmholtz

    module pure elemental subroutine calc_w_helmholtz(self, T_in, rho_in, w, prop_in)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: w
        type(type_iapws_phi_property), intent(inout), optional :: prop_in

        real(real64) :: tau, delta
        type(type_iapws_phi_property) :: props

        tau = self%T_c / T_in
        delta = rho_in / self%rho_c

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_phi(tau, delta, props)
        end if

        w = sqrt(self%R * T_in * &
                 (2.0d0 * delta * props%phi_d + delta**2 * props%phi_dd - &
                  ((delta * props%phi_d - delta * tau * props%phi_dt)**2) / &
                  (tau**2 * props%phi_tt)))

    end subroutine calc_w_helmholtz

end submodule iapws_property_helmholtz
