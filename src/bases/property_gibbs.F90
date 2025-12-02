submodule(module_iapws) property_gibbs
    implicit none

contains
    module pure elemental subroutine calc_properties_gibbs(self, T_in, p_in, property)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        type(type_iapws_property), intent(inout) :: property

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        call self%calc_gamma(tau, pi, props)

        call self%calc_nu(T_in, p_in, property%nu, props)
        call self%calc_rho(T_in, p_in, property%rho, props)
        call self%calc_u(T_in, p_in, property%u, props)
        call self%calc_h(T_in, p_in, property%h, props)
        call self%calc_s(T_in, p_in, property%s, props)
        call self%calc_cp(T_in, p_in, property%cp, props)
        call self%calc_cv(T_in, p_in, property%cv, props)
        call self%calc_w(T_in, p_in, property%w, props)

    end subroutine calc_properties_gibbs

    module pure elemental subroutine calc_nu_gibbs(self, T_in, p_in, nu, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: nu
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        nu = (pi * props%gamma_p) * self%R * T_in / p_in
    end subroutine calc_nu_gibbs

    module pure elemental subroutine calc_rho_gibbs(self, T_in, p_in, rho, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: rho
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props
        real(real64) :: nu

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        call self%calc_nu(T_in, p_in, nu, props)
        rho = 1.0d0 / nu
    end subroutine calc_rho_gibbs

    module pure elemental subroutine calc_u_gibbs(self, T_in, p_in, u, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: u
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        u = (tau * props%gamma_t - pi * props%gamma_p) * self%R * T_in

    end subroutine calc_u_gibbs

    module pure elemental subroutine calc_s_gibbs(self, T_in, p_in, s, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: s
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        s = (tau * props%gamma_t - props%gamma) * self%R

    end subroutine calc_s_gibbs

    module pure elemental subroutine calc_h_gibbs(self, T_in, p_in, h, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: h
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if
        h = tau * props%gamma_t * self%R * T_in

    end subroutine calc_h_gibbs

    module pure elemental subroutine calc_cp_gibbs(self, T_in, p_in, cp, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: cp
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        cp = (-tau**2 * props%gamma_tt) * self%R

    end subroutine calc_cp_gibbs

    module pure elemental subroutine calc_cv_gibbs(self, T_in, p_in, cv, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: cv
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        cv = (-tau**2 * props%gamma_tt - ((props%gamma_p - tau * props%gamma_pt)**2) / props%gamma_pp) * self%R
    end subroutine calc_cv_gibbs

    module pure elemental subroutine calc_w_gibbs(self, T_in, p_in, w, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: w
        type(type_iapws_gamma_property), intent(in), optional :: prop_in

        real(real64) :: tau, pi
        real(real64) :: w_sq
        type(type_iapws_gamma_property) :: props

        tau = self%T_star / T_in
        pi = p_in / self%p_star

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(tau, pi, props)
        end if

        w_sq = (props%gamma_p**2 / ((props%gamma_p - tau * props%gamma_pt)**2 / &
                                    (tau**2 * props%gamma_tt) - props%gamma_pp)) * self%R * T_in

        w = sqrt(max(w_sq, 0.0d0))

    end subroutine calc_w_gibbs

end submodule property_gibbs
