submodule(iapws) iapws_base
    implicit none
contains

    module pure elemental subroutine calc_properties_helmholtz(self, T_in, rho_in, property)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        type(type_iapws_property), intent(inout) :: property

    end subroutine calc_properties_helmholtz

    module pure elemental subroutine calc_pressure_helmholtz(self, T_in, rho_in, pressure)
        implicit none
        class(abst_iapws_helmholtz), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: rho_in
        real(real64), intent(inout) :: pressure

        type(type_iapws_phi_property) :: props

        call self%calc_phi(T_in, rho_in, props)

    end subroutine calc_pressure_helmholtz
end submodule iapws_base
