module batch_extension_model
    use iso_fortran_env,only:real64
    use iapws,only:type_iapws06,type_iapws97,type_iapws_property
    use module_iapws,only:type_iapws_gibbs_coefficient
    implicit none
    type,extends(type_iapws97)::custom_water
    contains
        procedure::calc_properties=>custom_properties
    end type
    type,extends(type_iapws06)::custom_ice
    contains
        procedure::calc_cp=>custom_cp
    end type
contains
    pure elemental subroutine custom_properties(self,T_in,p_in,property)
        class(custom_water),intent(in)::self
        real(real64),intent(in)::T_in,p_in
        type(type_iapws_property),intent(inout)::property
        call self%type_iapws97%calc_properties(T_in,p_in,property)
        property%cp=1777d0
    end subroutine
    pure elemental subroutine custom_cp(self,T_in,P_in,cp,prop_in)
        class(custom_ice),intent(in)::self
        real(real64),intent(in)::T_in,P_in
        real(real64),intent(inout)::cp
        type(type_iapws_gibbs_coefficient),intent(in),optional::prop_in
        cp=1234d0
    end subroutine
end module
program batch_extension
    use batch_extension_model
    implicit none
    type(custom_ice)::model
    type(custom_water)::water
    real(real64)::t(35),p(35),cp(35)
    t=260d0;p=1d5
    call model%initialize()
    call model%calc_properties_batch(t,p,cp=cp)
    if(any(abs(cp-1234d0)>1d-12)) error stop 'batch must respect overridden property methods'
    t=300d0
    call water%initialize()
    call water%calc_properties_batch(t,p,cp=cp)
    if(any(abs(cp-1777d0)>1d-12)) error stop 'batch must respect overridden IF97 aggregate' 
end program
