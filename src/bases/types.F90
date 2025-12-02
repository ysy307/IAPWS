submodule(module_iapws) iapws_types
    implicit none
contains
    module pure elemental subroutine reset_iapws_phi_property(self)
        implicit none
        !> IAPWS Helmholtz properties object
        class(type_iapws_phi_property), intent(inout) :: self

        self%phi = 0.0d0
        self%phi_d = 0.0d0
        self%phi_t = 0.0d0
        self%phi_dd = 0.0d0
        self%phi_dt = 0.0d0
        self%phi_tt = 0.0d0

        self%phi0 = 0.0d0
        self%phi0_d = 0.0d0
        self%phi0_t = 0.0d0
        self%phi0_dd = 0.0d0
        self%phi0_dt = 0.0d0
        self%phi0_tt = 0.0d0
        self%phir = 0.0d0
        self%phir_d = 0.0d0
        self%phir_t = 0.0d0
        self%phir_dd = 0.0d0
        self%phir_dt = 0.0d0
        self%phir_tt = 0.0d0
    end subroutine reset_iapws_phi_property

    module pure elemental subroutine reset_iapws_gamma_property(self)
        implicit none
        !> IAPWS Gibbs properties object
        class(type_iapws_gamma_property), intent(inout) :: self

        self%gamma = 0.0d0
        self%gamma_p = 0.0d0
        self%gamma_t = 0.0d0
        self%gamma_pp = 0.0d0
        self%gamma_tt = 0.0d0
        self%gamma_pt = 0.0d0

        self%gamma0 = 0.0d0
        self%gamma0_p = 0.0d0
        self%gamma0_t = 0.0d0
        self%gamma0_pp = 0.0d0
        self%gamma0_tt = 0.0d0
        self%gamma0_pt = 0.0d0
        self%gammar = 0.0d0
        self%gammar_p = 0.0d0
        self%gammar_t = 0.0d0
        self%gammar_pp = 0.0d0
        self%gammar_tt = 0.0d0
        self%gammar_pt = 0.0d0
    end subroutine reset_iapws_gamma_property

end submodule iapws_types
