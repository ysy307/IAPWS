submodule(module_iapws) property_gibbs
    implicit none

contains
    module pure elemental subroutine calc_properties_gibbs(self, T_in, P_in, property)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        type(type_iapws_property), intent(inout) :: property

        type(type_iapws_gibbs_coefficient) :: coef

        call self%calc_gamma(T_in, P_in, coef)

        call self%calc_nu(T_in, P_in, property%nu, coef)
        call self%calc_rho(T_in, P_in, property%rho, coef)
        call self%calc_u(T_in, P_in, property%u, coef)
        call self%calc_h(T_in, P_in, property%h, coef)
        call self%calc_s(T_in, P_in, property%s, coef)
        call self%calc_cp(T_in, P_in, property%cp, coef)
        call self%calc_cv(T_in, P_in, property%cv, coef)
        call self%calc_w(T_in, P_in, property%w, coef)
        call self%calc_alpha(T_in, P_in, property%alpha, coef)
        call self%calc_beta(T_in, P_in, property%beta, coef)
        call self%calc_kappa_T(T_in, P_in, property%kappa_T, coef)
        call self%calc_kappa_s(T_in, P_in, property%kappa_s, coef)

    end subroutine calc_properties_gibbs

    module pure elemental subroutine calc_nu_gibbs(self, T_in, P_in, nu, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: nu
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        nu = props%g_p
    end subroutine calc_nu_gibbs

    module pure elemental subroutine calc_rho_gibbs(self, T_in, P_in, rho, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: rho
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props
        real(real64) :: nu

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        call self%calc_nu(T_in, P_in, nu, props)
        rho = 1.0d0 / nu
    end subroutine calc_rho_gibbs

    module pure elemental subroutine calc_u_gibbs(self, T_in, P_in, u, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: u
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        u = props%g - T_in * props%g_t - P_in * props%g_p

    end subroutine calc_u_gibbs

    module pure elemental subroutine calc_s_gibbs(self, T_in, P_in, s, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: s
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        s = -props%g_t

    end subroutine calc_s_gibbs

    module pure elemental subroutine calc_h_gibbs(self, T_in, P_in, h, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: h
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        h = props%g - T_in * props%g_t

    end subroutine calc_h_gibbs

    module pure elemental subroutine calc_cp_gibbs(self, T_in, P_in, cp, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: cp
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        cp = -T_in * props%g_tt

    end subroutine calc_cp_gibbs

    module pure elemental subroutine calc_cv_gibbs(self, T_in, P_in, cv, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: cv
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        cv = -T_in * props%g_tt + T_in * (props%g_pt**2) / props%g_pp
    end subroutine calc_cv_gibbs

    module pure elemental subroutine calc_w_gibbs(self, T_in, P_in, w, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: w
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        real(real64) :: w_sq
        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        w_sq = (props%g_p**2) / ((props%g_pt**2 / props%g_tt) - props%g_pp)
        w = sqrt(max(w_sq, 0.0d0))

    end subroutine calc_w_gibbs

    module pure elemental subroutine calc_alpha_gibbs(self, T_in, P_in, alpha, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: alpha
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in
        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        alpha = props%g_pt / props%g_p

    end subroutine calc_alpha_gibbs

    module pure elemental subroutine calc_beta_gibbs(self, T_in, P_in, beta, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: beta
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in

        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        beta = -props%g_pt / props%g_pp

    end subroutine calc_beta_gibbs

    module pure elemental subroutine calc_kappa_T_gibbs(self, T_in, P_in, kappa_T, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: kappa_T
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in
        type(type_iapws_gibbs_coefficient) :: props

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        kappa_T = -props%g_pp / props%g_p

    end subroutine calc_kappa_T_gibbs

    module pure elemental subroutine calc_kappa_s_gibbs(self, T_in, P_in, kappa_s, prop_in)
        implicit none
        class(abst_iapws_gibbs), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: P_in
        real(real64), intent(inout) :: kappa_s
        type(type_iapws_gibbs_coefficient), intent(in), optional :: prop_in
        type(type_iapws_gibbs_coefficient) :: props
        real(real64) :: term

        if (present(prop_in)) then
            props = prop_in
        else
            call self%calc_gamma(T_in, P_in, props)
        end if

        kappa_s = (props%g_pt**2 - props%g_tt * props%g_pp) / (props%g_p * props%g_tt)

    end subroutine calc_kappa_s_gibbs

end submodule property_gibbs
