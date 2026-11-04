program saturation_derivatives
    use iso_fortran_env, only: real64
    use ieee_arithmetic
    use iapws, only: type_iapws97
    implicit none
    type(type_iapws97) :: water
    real(real64), parameter :: ts(8)=[273.15d0,300d0,500d0,623.15d0,623.16d0,630d0,640d0,647d0]
    real(real64) :: dc_v(8),dc_l(8),dr_v(8),dr_l(8),dh(8),rv,rl,cv,cl,vp,lp,vm,lm,hp,hm,h
    real(real64) :: bad(6),pv(6),pl(6),dt,dp,rp,rm,actual
    real(real64) :: bulk(8,5)
    integer :: i,j
    call water%initialize()
    call water%calc_saturation_dcp_dT(ts,dc_v,dc_l)
    call water%calc_saturation_density_dT(ts,dr_v,dr_l)
    call water%calc_latent_heat_dT(ts,dh)
    call water%calc_saturation_properties(ts,dcp_dT_vap=bulk(:,1),dcp_dT_liq=bulk(:,2), &
        density_dT_vap=bulk(:,3),density_dT_liq=bulk(:,4),dlatent_heat_dT=bulk(:,5))
    do i=1,size(ts)
        call check(bulk(i,1),dc_v(i),1d-12,'bulk vapor cp compatibility')
        call check(bulk(i,2),dc_l(i),1d-12,'bulk liquid cp compatibility')
        call check(bulk(i,3),dr_v(i),1d-12,'bulk vapor density compatibility')
        call check(bulk(i,4),dr_l(i),1d-12,'bulk liquid density compatibility')
        call check(bulk(i,5),dh(i),1d-12,'bulk latent heat compatibility')
    end do
    dc_v=bulk(:,1); dc_l=bulk(:,2); dr_v=bulk(:,3); dr_l=bulk(:,4); dh=bulk(:,5)
    do i=1,size(ts)
        call water%calc_saturation_density(ts(i),rv,rl)
        if (.not.ieee_is_finite(rv).or..not.ieee_is_finite(rl)) error stop 'nonfinite phase density'
        if (rv>=rl) error stop 'distinct vapor/liquid roots required'
        do j=1,2
            h=0.0002d0/real(j,real64)
            ! Stay on the selected branch at the lower endpoint and region join.
            if (i==1.or.i==4) then
                if (i==1) h=-h
                call water%calc_saturation_cp(ts(i),cv,cl)
                call water%calc_saturation_cp(ts(i)-h,vm,lm)
                vp=cv; lp=cl
                dt=h
                call water%calc_saturation_density(ts(i),rp,rm)
                call water%calc_saturation_density(ts(i)-h,rv,rl)
                call check(dr_v(i),(rp-rv)/h,1d-3,'density vapor join')
                call check(dr_l(i),(rm-rl)/h,1d-3,'density liquid join')
                call water%calc_latent_heat(hp,T_in=ts(i))
                call water%calc_latent_heat(hm,T_in=ts(i)-h)
            else
                dt=2*h
                call water%calc_saturation_cp(ts(i)+h,vp,lp)
                call water%calc_saturation_cp(ts(i)-h,vm,lm)
                call water%calc_saturation_density(ts(i)+h,rp,rm)
                call water%calc_saturation_density(ts(i)-h,rv,rl)
                call check(dr_v(i),(rp-rv)/dt,1d-3,'density vapor')
                call check(dr_l(i),(rm-rl)/dt,1d-3,'density liquid')
                call water%calc_latent_heat(hp,T_in=ts(i)+h)
                call water%calc_latent_heat(hm,T_in=ts(i)-h)
            end if
            call check(dc_v(i),(vp-vm)/dt,1d-3,'cp vapor')
            call check(dc_l(i),(lp-lm)/dt,1d-3,'cp liquid')
            call check(dh(i),(hp-hm)/dt,1d-3,'latent heat')
        end do
        call water%calc_saturation_dcp_dT(ts(i),dcp_dT_liq=actual)
        call check(actual,dc_l(i),1d-12,'optional liquid cp')
        call water%calc_saturation_density_dT(ts(i),drho_dT_vap=actual)
        call check(actual,dr_v(i),1d-12,'optional vapor density')
    end do
    bad=[273.14d0,647.096d0,648d0,-1d0,ieee_value(0d0,ieee_quiet_nan),ieee_value(0d0,ieee_positive_inf)]
    call water%calc_saturation_dcp_dT(bad,pv,pl)
    if (.not.all(ieee_is_nan(pv)).or..not.all(ieee_is_nan(pl))) error stop 'invalid saturation cp'
    call water%calc_saturation_density_dT(bad,pv,pl)
    if (.not.all(ieee_is_nan(pv)).or..not.all(ieee_is_nan(pl))) error stop 'invalid saturation density'
    call water%calc_latent_heat_dT(bad,pv)
    if (.not.all(ieee_is_nan(pv))) error stop 'invalid latent derivative'
    call water%calc_saturation_density(647.096d0,rv,rl)
    call check(rv,322d0,1d-12,'critical vapor density')
    call check(rl,322d0,1d-12,'critical liquid density')
    call water%calc_latent_heat(actual,T_in=647.096d0)
    if (actual/=0d0) error stop 'critical latent heat'
    ! Independent EOS/root investigation: fits cease to bracket a vapor root
    ! extremely close to Tc. A failed phase pair must not expose a partial state.
    call water%calc_saturation_density(647.09599d0,rv,rl)
    if (.not.ieee_is_nan(rv).or..not.ieee_is_nan(rl)) error stop 'unbracketed phase pair must be NaN'
    ! Fixed-pressure/temperature density partial regression in region 3.
    h=0.001d0; dp=100d0
    call water%calc_drho_dT(650d0,25d6,actual)
    call water%calc_rho(650d0+h,25d6,rp)
    call water%calc_rho(650d0-h,25d6,rm)
    call check(actual,(rp-rm)/(2*h),1d-5,'region3 drho/dT')
    call water%calc_drho_dp(650d0,25d6,actual)
    call water%calc_rho(650d0,25d6+dp,rp)
    call water%calc_rho(650d0,25d6-dp,rm)
    call check(actual,(rp-rm)/(2*dp),1d-5,'region3 drho/dp')
contains
    subroutine check(actual,expected,tol,label)
        real(real64),intent(in)::actual,expected,tol
        character(*),intent(in)::label
        if (.not.ieee_is_finite(actual).or..not.ieee_is_finite(expected)) then
            print *,label,actual,expected
            error stop 'nonfinite saturation result'
        end if
        if(abs(actual-expected)>tol*max(abs(expected),1d-7)) then
            print *,label,actual,expected
            error stop 'saturation derivative mismatch'
        end if
    end subroutine
end program
