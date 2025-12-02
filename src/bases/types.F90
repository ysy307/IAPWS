submodule(module_iapws) iapws_types
    implicit none
contains
    module pure elemental subroutine reset_iapws_helmholtz_property(self)
        implicit none
        !> IAPWS Helmholtz properties object
        class(type_iapws_helmholtz_property), intent(inout) :: self

        self%f = 0.0d0
        self%f_d = 0.0d0
        self%f_t = 0.0d0
        self%f_dd = 0.0d0
        self%f_dt = 0.0d0
        self%f_tt = 0.0d0

        self%f0 = 0.0d0
        self%f0_d = 0.0d0
        self%f0_t = 0.0d0
        self%f0_dd = 0.0d0
        self%f0_dt = 0.0d0
        self%f0_tt = 0.0d0
        self%fr = 0.0d0
        self%fr_d = 0.0d0
        self%fr_t = 0.0d0
        self%fr_dd = 0.0d0
        self%fr_dt = 0.0d0
        self%fr_tt = 0.0d0
    end subroutine reset_iapws_helmholtz_property

    module pure elemental subroutine reset_iapws_gibbs_coefficient(self)
        implicit none
        !> IAPWS Gibbs properties object
        class(type_iapws_gibbs_coefficient), intent(inout) :: self

        self%g = 0.0d0
        self%g_p = 0.0d0
        self%g_t = 0.0d0
        self%g_pp = 0.0d0
        self%g_tt = 0.0d0
        self%g_pt = 0.0d0

        self%g0 = 0.0d0
        self%g0_p = 0.0d0
        self%g0_t = 0.0d0
        self%g0_pp = 0.0d0
        self%g0_tt = 0.0d0
        self%g0_pt = 0.0d0
        self%gr = 0.0d0
        self%gr_p = 0.0d0
        self%gr_t = 0.0d0
        self%gr_pp = 0.0d0
        self%gr_tt = 0.0d0
        self%gr_pt = 0.0d0
    end subroutine reset_iapws_gibbs_coefficient
end submodule iapws_types
