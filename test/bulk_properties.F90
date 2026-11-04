program bulk_properties
    use iso_fortran_env, only: real64, int64
    use ieee_arithmetic
    use iapws, only: type_iapws97
    implicit none
    type(type_iapws97) :: w
    integer, parameter :: n=271
    real(real64) :: ts(n), a(n,14), b(n,14), single(n)
    real(real64) :: tc(24),pc(24),ct(24),cp(24),ref_ct(24),ref_cp(24),nan,inf
    character(24), parameter :: labels(14)=[character(24) :: &
        'rho_vap','rho_liq','cp_vap','cp_liq','drho_dT_vap','drho_dT_liq', &
        'drho_dP_vap','drho_dP_liq','density_dT_vap','density_dT_liq', &
        'dcp_dT_vap','dcp_dT_liq','latent_heat','dlatent_heat_dT']
    integer :: i,j
    call w%initialize()
    nan=ieee_value(0d0,ieee_quiet_nan)
    inf=ieee_value(0d0,ieee_positive_inf)
    ts(:16)=[273.15d0,300d0,500d0,623.15d0,623.16d0,640d0,647d0,647.096d0, &
        270d0,648d0,nan,inf,-inf,647.09599d0, &
        nearest(273.15d0,-1d0),nearest(647.096d0,1d0)]
    do i=17,n
        ts(i)=273.15d0+(647.09d0-273.15d0)*real(i-17,real64)/real(n-17,real64)
    end do
    call w%calc_saturation_properties(ts,rho_vap=a(:,1),rho_liq=a(:,2), &
        cp_vap=a(:,3),cp_liq=a(:,4),drho_dT_vap=a(:,5),drho_dT_liq=a(:,6), &
        drho_dP_vap=a(:,7),drho_dP_liq=a(:,8),density_dT_vap=a(:,9),density_dT_liq=a(:,10), &
        dcp_dT_vap=a(:,11),dcp_dT_liq=a(:,12),latent_heat=a(:,13),dlatent_heat_dT=a(:,14))
    call w%calc_saturation_density(ts,b(:,1),b(:,2))
    call w%calc_saturation_cp(ts,b(:,3),b(:,4))
    call w%calc_saturation_drho_dT(ts,b(:,5),b(:,6))
    call w%calc_saturation_drho_dP(ts,b(:,7),b(:,8))
    call w%calc_saturation_density_dT(ts,b(:,9),b(:,10))
    call w%calc_saturation_dcp_dT(ts,b(:,11),b(:,12))
    call w%calc_latent_heat(b(:,13),T_in=ts)
    call w%calc_latent_heat_dT(ts,b(:,14))
    do j=1,14
        do i=1,size(ts)
            call same(a(i,j),b(i,j),labels(j),ts(i))
        end do
    end do
    if(any(.not.ieee_is_finite(a(:7,:)))) error stop 'finite bulk saturation states'
    if(any(.not.ieee_is_finite(a(17:,:)))) error stop 'finite dense bulk saturation grid'
    ! Rows 9:16 are invalid or the documented failed root at 647.09599 K.
    if(any(.not.ieee_is_nan(a(9:16,:)))) error stop 'invalid saturation states must return NaN'
    call same(a(8,1),322d0,'critical vapor density',ts(8))
    call same(a(8,2),322d0,'critical liquid density',ts(8))
    call same(a(8,13),0d0,'critical latent heat',ts(8))
    if(any(.not.ieee_is_nan(a(8,3:12))).or..not.ieee_is_nan(a(8,14))) error stop 'critical bulk singularities'
    ! Each output must work independently: in particular, value-only and
    ! derivative-only requests must preserve critical/invalid-state semantics.
    call w%calc_saturation_properties(ts,rho_vap=single)
    do i=1,n
        call same(single(i),b(i,1),labels(1),ts(i))
    end do
    call w%calc_saturation_properties(ts,rho_liq=single)
    do i=1,n
        call same(single(i),b(i,2),labels(2),ts(i))
    end do
    call w%calc_saturation_properties(ts,cp_vap=single)
    do i=1,n
        call same(single(i),b(i,3),labels(3),ts(i))
    end do
    call w%calc_saturation_properties(ts,cp_liq=single)
    do i=1,n
        call same(single(i),b(i,4),labels(4),ts(i))
    end do
    call w%calc_saturation_properties(ts,drho_dT_vap=single)
    do i=1,n
        call same(single(i),b(i,5),labels(5),ts(i))
    end do
    call w%calc_saturation_properties(ts,drho_dT_liq=single)
    do i=1,n
        call same(single(i),b(i,6),labels(6),ts(i))
    end do
    call w%calc_saturation_properties(ts,drho_dP_vap=single)
    do i=1,n
        call same(single(i),b(i,7),labels(7),ts(i))
    end do
    call w%calc_saturation_properties(ts,drho_dP_liq=single)
    do i=1,n
        call same(single(i),b(i,8),labels(8),ts(i))
    end do
    call w%calc_saturation_properties(ts,density_dT_vap=single)
    do i=1,n
        call same(single(i),b(i,9),labels(9),ts(i))
    end do
    call w%calc_saturation_properties(ts,density_dT_liq=single)
    do i=1,n
        call same(single(i),b(i,10),labels(10),ts(i))
    end do
    call w%calc_saturation_properties(ts,dcp_dT_vap=single)
    do i=1,n
        call same(single(i),b(i,11),labels(11),ts(i))
    end do
    call w%calc_saturation_properties(ts,dcp_dT_liq=single)
    do i=1,n
        call same(single(i),b(i,12),labels(12),ts(i))
    end do
    call w%calc_saturation_properties(ts,latent_heat=single)
    do i=1,n
        call same(single(i),b(i,13),labels(13),ts(i))
    end do
    call w%calc_saturation_properties(ts,dlatent_heat_dT=single)
    do i=1,n
        call same(single(i),b(i,14),labels(14),ts(i))
    end do
    ! A value-only batch, a liquid-only batch, and the six-output client case.
    call w%calc_saturation_properties(ts,rho_vap=a(:,1),rho_liq=a(:,2), &
        cp_vap=a(:,3),cp_liq=a(:,4),latent_heat=a(:,13))
    call same_columns([1,2,3,4,13])
    call w%calc_saturation_properties(ts,rho_liq=a(:,2),cp_liq=a(:,4), &
        density_dT_liq=a(:,10),dcp_dT_liq=a(:,12))
    call same_columns([2,4,10,12])
    call w%calc_saturation_properties(ts,rho_vap=a(:,1),density_dT_vap=a(:,9), &
        cp_vap=a(:,3),dcp_dT_vap=a(:,11),latent_heat=a(:,13),dlatent_heat_dT=a(:,14))
    do j=1,14
        do i=1,n
            call same(a(i,j),b(i,j),labels(j),ts(i))
        end do
    end do
    ! Pure elemental bulk dispatch: valid/invalid elements must be independent.
    tc=[300d0,700d0,650d0,1500d0,243.15d0,1073.15d0,2273.15d0, &
        623.15d0,623.16d0,640d0,0d0,-1d0,243.14d0,2274d0, &
        300d0,300d0,300d0,1500d0,nan,inf,-inf,300d0,300d0,647.096d0]
    pc=[3d6,3d6,25d6,30d6,1d5,100d6,50d6,20d6,20d6,20d6, &
        1d5,1d5,1d5,1d5,0d0,-1d0,101d6,51d6,1d5,1d5,1d5,nan,inf,22.064d6]
    call evaluate_pure(w,tc,pc,ct,cp)
    call w%calc_dcp_dT(tc,pc,ref_ct)
    call w%calc_dcp_dp(tc,pc,ref_cp)
    do i=1,size(tc)
        call same(ct(i),ref_ct(i),'ordinary cp temperature',tc(i))
        call same(cp(i),ref_cp(i),'ordinary cp pressure',tc(i))
        if(i<=10) then
            if(.not.ieee_is_finite(ct(i)).or..not.ieee_is_finite(cp(i))) error stop 'valid bulk cp partials'
        else
            if(.not.ieee_is_nan(ct(i)).or..not.ieee_is_nan(cp(i))) error stop 'invalid bulk cp partials'
        end if
    end do
    ! Scalar temperature broadcast against an array of pressures.
    call w%calc_cp_derivatives(300d0,pc,ct,cp)
    call w%calc_dcp_dT(300d0,pc,ref_ct)
    call w%calc_dcp_dp(300d0,pc,ref_cp)
    do i=1,size(pc)
        call same(ct(i),ref_ct(i),'broadcast cp temperature',300d0)
        call same(cp(i),ref_cp(i),'broadcast cp pressure',300d0)
    end do
#ifdef _OPENMP
    call check_parallel_properties()
#endif
contains
    subroutine same_columns(columns)
        integer,intent(in)::columns(:)
        integer::row,col,index
        do index=1,size(columns)
            col=columns(index)
            do row=1,n
                call same(a(row,col),b(row,col),labels(col),ts(row))
            end do
        end do
    end subroutine

    pure subroutine evaluate_pure(model,t,p,ct,cp)
        type(type_iapws97),intent(in)::model
        real(real64),intent(in)::t(:),p(:)
        real(real64),intent(out)::ct(:),cp(:)
        call model%calc_cp_derivatives(t,p,ct,cp)
    end subroutine

    subroutine same(x,y,label,t)
        real(real64),intent(in)::x,y
        character(*),intent(in)::label
        real(real64),intent(in)::t
        if(ieee_is_nan(x).and.ieee_is_nan(y)) return
        if(transfer(x,0_int64)/=transfer(y,0_int64)) then
            print *, label,t,x,y
            error stop 'bulk/scalar mismatch'
        end if
    end subroutine
end program

#ifdef _OPENMP
! One immutable model can be shared by independent workers. Check pointwise
! results against serial evaluation, with reversed scheduling and mixed regions.
subroutine check_parallel_properties
    use iso_fortran_env, only: real64, int64
    use omp_lib, only: omp_get_num_threads
    use ieee_arithmetic, only: ieee_is_finite
    use iapws, only: type_iapws97, type_iapws_property
    implicit none
    integer,parameter :: n=256, outputs=25
    type(type_iapws97) :: water
    real(real64) :: expected(outputs,n),actual(outputs,n)
    integer :: i,threads
    call water%initialize()
    do i=1,n
        call evaluate(i,expected(:,i))
    end do
    threads=0
    !$omp parallel default(none) shared(water,actual,threads) private(i) num_threads(4)
    !$omp single
    threads=omp_get_num_threads()
    !$omp end single
    !$omp do schedule(dynamic,1)
    do i=n,1,-1
        call evaluate(i,actual(:,i))
    end do
    !$omp end do
    !$omp end parallel
    if(threads<2) error stop 'parallel test needs multiple workers'
    if(any(transfer(expected,[0_int64])/=transfer(actual,[0_int64]))) error stop 'parallel property mismatch'
    print *, 'Mixed-region properties match serial results with workers:',threads
contains
    subroutine evaluate(index,v)
        integer,intent(in)::index
        real(real64),intent(out)::v(outputs)
        type(type_iapws_property)::props
        real(real64)::t,p,ts
        select case(mod(index,4))
        case(0)
            t=280d0+mod(index,100)*3d0; p=30d6
        case(1)
            t=700d0+mod(index,100)*3d0; p=1d5
        case(2)
            t=650d0+mod(index,100)*0.4d0; p=30d6
        case(3)
            t=1200d0+mod(index,100)*9d0; p=10d6
        end select
        call water%calc_properties(t,p,props)
        v(:9)=[props%rho,props%cp,props%cv,props%h,props%u,props%s,props%w,props%drho_dT,props%drho_dP]
        call water%calc_cp_derivatives(t,p,v(10),v(11))
        ts=273.15d0+mod(index,256)*1.46d0
        call water%calc_saturation_properties(ts,rho_vap=v(12),rho_liq=v(13),cp_vap=v(14),cp_liq=v(15), &
            drho_dT_vap=v(16),drho_dT_liq=v(17),drho_dP_vap=v(18),drho_dP_liq=v(19), &
            density_dT_vap=v(20),density_dT_liq=v(21),dcp_dT_vap=v(22),dcp_dT_liq=v(23), &
            latent_heat=v(24),dlatent_heat_dT=v(25))
        if(.not.all(ieee_is_finite(v))) error stop 'nonfinite parallel property'
    end subroutine
end subroutine check_parallel_properties
#endif
