submodule(iapws) iapws_base
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

end submodule iapws_base
