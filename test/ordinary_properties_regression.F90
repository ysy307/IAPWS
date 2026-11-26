! Cross-check the public ordinary-property path against the already independent
! saturation fixtures and consistency relations. Temperature/pressure points
! include both subcritical branches, supercritical states and the region join.
program ordinary_properties_regression
    use iso_fortran_env, only: real64, int64
    use ieee_arithmetic, only: ieee_is_finite
    use iapws, only: type_iapws97, type_iapws_property
    implicit none
    type(type_iapws97) :: water
    type(type_iapws_property) :: props
    real(real64),parameter :: temperatures(12)=[623.16d0,627.15d0,631.15d0,635.15d0, &
        639.15d0,640d0,646.15d0,650d0,650d0,675d0,700d0,750d0]
    real(real64),parameter :: pressures(12)=[20d6,17d6,18d6,19d6,20d6,20d6, &
        22d6,25d6,30d6,30d6,30d6,80d6]
    ! Densities from the bracketed forward solve at 9b0dfef.
    real(real64), parameter :: density_ref(12) = [ &
        600.6170143181862d0, 114.93549968338114d0, 129.41312747278178d0, 146.60623318331273d0, &
        167.90922446153317d0, 160.57788700157562d0, 433.57024154875256d0, 488.87505207909874d0, &
        549.58058941471359d0, 332.81815541638588d0, 184.18016875974061d0, 506.66424785058132d0]
    real(real64)::rho,cp,rho_T,rho_p,ct,cp_p,ht,hp,plus,minus
    integer::i
    call water%initialize()
    do i=1,size(temperatures)
        call water%calc_properties(temperatures(i),pressures(i),props)
        call water%calc_rho(temperatures(i),pressures(i),rho)
        call water%calc_cp(temperatures(i),pressures(i),cp)
        call same(rho,props%rho)
        call close(rho,density_ref(i),2048d0*epsilon(rho))
        call close(props%p,pressures(i),2048d0*epsilon(rho))
        call same(cp,props%cp)
        call water%calc_drho_dT(temperatures(i),pressures(i),rho_T)
        call water%calc_drho_dp(temperatures(i),pressures(i),rho_p)
        call close(rho_T,props%drho_dT,1d-12)
        call close(rho_p,props%drho_dP,1d-12)
        call water%calc_cp_derivatives(temperatures(i),pressures(i),ct,cp_p)
        ht=0.001d0; hp=50d0
        call water%calc_cp(temperatures(i)+ht,pressures(i),plus)
        call water%calc_cp(temperatures(i)-ht,pressures(i),minus)
        call close(ct,(plus-minus)/(2*ht),1d-4)
        call water%calc_cp(temperatures(i),pressures(i)+hp,plus)
        call water%calc_cp(temperatures(i),pressures(i)-hp,minus)
        call close(cp_p,(plus-minus)/(2*hp),1d-4)
    end do
contains
    subroutine same(a,b)
        real(real64),intent(in)::a,b
        if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) error stop 'nonfinite ordinary properties'
        if(transfer(a,0_int64)/=transfer(b,0_int64)) error stop 'ordinary property consistency'
    end subroutine
    subroutine close(a,b,tolerance)
        real(real64),intent(in)::a,b,tolerance
        if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) error stop 'nonfinite ordinary derivative'
        if(abs(a-b)>tolerance*max(abs(a),abs(b),1d-10)) error stop 'ordinary derivative consistency'
    end subroutine
end program
