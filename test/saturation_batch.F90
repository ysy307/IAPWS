program saturation_batch
    use iso_fortran_env,only:real64
    use ieee_arithmetic
    use iapws,only:type_iapws97
    implicit none
    integer,parameter::n=137
    type(type_iapws97)::water
    real(real64)::t(n),v(n,14),ref(14),single(n),ts(2*n),out(2*n)
    integer::i,j,status(n)
    call water%initialize()
    do i=1,n
        t(i)=273.15d0+350d0*real(mod(i*97,137),real64)/137d0
    end do
    t(:8)=[273.15d0,623.15d0,623.1501d0,640d0,647.096d0,647.09599d0, &
        270d0,ieee_value(0d0,ieee_quiet_nan)]
    call water%calc_saturation_properties_batch(t,rho_vap=v(:,1),rho_liq=v(:,2), &
        cp_vap=v(:,3),cp_liq=v(:,4),drho_dT_vap=v(:,5),drho_dT_liq=v(:,6), &
        drho_dP_vap=v(:,7),drho_dP_liq=v(:,8),density_dT_vap=v(:,9),density_dT_liq=v(:,10), &
        dcp_dT_vap=v(:,11),dcp_dT_liq=v(:,12),latent_heat=v(:,13),dlatent_heat_dT=v(:,14),status=status)
    do i=1,n
        call water%calc_saturation_properties(t(i),rho_vap=ref(1),rho_liq=ref(2), &
            cp_vap=ref(3),cp_liq=ref(4),drho_dT_vap=ref(5),drho_dT_liq=ref(6), &
            drho_dP_vap=ref(7),drho_dP_liq=ref(8),density_dT_vap=ref(9),density_dT_liq=ref(10), &
            dcp_dT_vap=ref(11),dcp_dT_liq=ref(12),latent_heat=ref(13),dlatent_heat_dT=ref(14))
        do j=1,14
            call check(v(i,j),ref(j))
        end do
    end do
    if(any(status(7:8)/=-1).or.any(status(:6)/=0)) error stop 'saturation batch status'
    call water%calc_saturation_properties_batch(t,rho_vap=single)
    do i=1,n
        call check(single(i),v(i,1))
    end do
    call water%calc_saturation_properties_batch(t,rho_liq=single)
    do i=1,n
        call check(single(i),v(i,2))
    end do
    call water%calc_saturation_properties_batch(t,cp_vap=single)
    do i=1,n
        call check(single(i),v(i,3))
    end do
    call water%calc_saturation_properties_batch(t,cp_liq=single)
    do i=1,n
        call check(single(i),v(i,4))
    end do
    call water%calc_saturation_properties_batch(t,drho_dT_vap=single)
    do i=1,n
        call check(single(i),v(i,5))
    end do
    call water%calc_saturation_properties_batch(t,drho_dT_liq=single)
    do i=1,n
        call check(single(i),v(i,6))
    end do
    call water%calc_saturation_properties_batch(t,drho_dP_vap=single)
    do i=1,n
        call check(single(i),v(i,7))
    end do
    call water%calc_saturation_properties_batch(t,drho_dP_liq=single)
    do i=1,n
        call check(single(i),v(i,8))
    end do
    call water%calc_saturation_properties_batch(t,density_dT_vap=single)
    do i=1,n
        call check(single(i),v(i,9))
    end do
    call water%calc_saturation_properties_batch(t,density_dT_liq=single)
    do i=1,n
        call check(single(i),v(i,10))
    end do
    call water%calc_saturation_properties_batch(t,dcp_dT_vap=single)
    do i=1,n
        call check(single(i),v(i,11))
    end do
    call water%calc_saturation_properties_batch(t,dcp_dT_liq=single)
    do i=1,n
        call check(single(i),v(i,12))
    end do
    call water%calc_saturation_properties_batch(t,latent_heat=single)
    do i=1,n
        call check(single(i),v(i,13))
    end do
    call water%calc_saturation_properties_batch(t,dlatent_heat_dT=single)
    do i=1,n
        call check(single(i),v(i,14))
    end do
    ts=ieee_value(0d0,ieee_quiet_nan);out=-999d0
    ts(1:2*n:2)=t
    call water%calc_saturation_properties_batch(ts(1:2*n:2),rho_vap=out(1:2*n:2))
    do i=1,n
        call check(out(2*i-1),v(i,1))
    end do
    if(any(out(2:2*n:2)/=-999d0)) error stop 'saturation output stride'
    t(:3)=[ieee_value(0d0,ieee_positive_inf),648d0,273.14d0]
    call water%calc_saturation_properties_batch(t,status=status)
    if(any(status(:3)/=-1)) error stop 'saturation metadata-only range failures'
    call water%calc_saturation_properties_batch(t(:0),rho_vap=single(:0),status=status(:0))
contains
    subroutine check(a,b)
        real(real64),intent(in)::a,b
        if(ieee_is_nan(a).and.ieee_is_nan(b)) return
        if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b)) error stop 'saturation batch nonfinite mismatch'
        if(abs(a-b)>5d-12*max(abs(a),abs(b),1d-10)) then
            print *,i,j,a,b
            error stop 'saturation batch scalar mismatch'
        end if
    end subroutine
end program
