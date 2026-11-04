program batch_shapes
    use iso_fortran_env,only:real64
    use iapws,only:type_iapws97,type_iapws95,type_iapws06
    implicit none
    type(type_iapws97)::water
    type(type_iapws95)::fluid
    type(type_iapws06)::ice
    real(real64)::t(3),p(3),cp(3)
    integer::status(3)
    character(24)::model,mode
    call get_command_argument(1,model);call get_command_argument(2,mode)
    call water%initialize();call fluid%initialize();call ice%initialize()
    t=300d0;p=1d5
    select case(trim(model))
    case('97')
        select case(trim(mode))
        case('input')
            call water%calc_properties_batch(t,p(:2),cp=cp)
        case('output')
            call water%calc_properties_batch(t,p,cp=cp(:2))
        case('status')
            call water%calc_properties_batch(t,p,status=status(:2))
        end select
    case('95')
        p=996d0
        select case(trim(mode))
        case('input')
            call fluid%calc_properties_batch(t,p(:2),cp=cp)
        case('output')
            call fluid%calc_properties_batch(t,p,cp=cp(:2))
        case('status')
            call fluid%calc_properties_batch(t,p,status=status(:2))
        end select
    case('sat')
        select case(trim(mode))
        case('output')
            call water%calc_saturation_properties_batch(t,rho_vap=cp(:2))
        case('status')
            call water%calc_saturation_properties_batch(t,status=status(:2))
        end select
    case('ice')
        t=260d0
        select case(trim(mode))
        case('input')
            call ice%calc_properties_batch(t,p(:2),cp=cp)
        case('output')
            call ice%calc_properties_batch(t,p,cp=cp(:2))
        case('status')
            call ice%calc_properties_batch(t,p,status=status(:2))
        end select
    end select
end program
