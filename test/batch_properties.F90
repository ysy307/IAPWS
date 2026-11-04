program batch_properties
    use iso_fortran_env,only:real64
    use ieee_arithmetic
    use iapws,only:type_iapws97,type_iapws95,type_iapws06,type_iapws_property
    implicit none
    integer,parameter::n=131
    type(type_iapws97)::water
    type(type_iapws95)::fluid
    type(type_iapws06)::ice
    type(type_iapws_property)::ref
    real(real64)::t(n),p(n),rho(n),cp(n),dt(n),dp(n),cv(n),w(n),u(n),s(n),h(n),nu(n),press(n)
    real(real64)::scratch(2*n),t_store(2*n),p_store(2*n),nan,inf
    integer::meta_store(2*n),region_store(2*n)
    integer::status(n),regions(n),i,k,m
    integer,parameter::sizes(12)=[0,1,3,4,7,8,31,32,33,63,64,n]
    call water%initialize();call fluid%initialize();call ice%initialize()
    do i=1,n
        select case(mod(i,5))
        case(0)
            t(i)=280d0+mod(i,100)*3d0;p(i)=30d6
        case(1)
            t(i)=700d0+mod(i,100)*3d0;p(i)=1d5
        case(2)
            t(i)=650d0+mod(i,100)*0.4d0;p(i)=30d6
        case(3)
            t(i)=1200d0+mod(i,100)*8d0;p(i)=10d6
        case(4)
            t(i)=640d0;p(i)=20d6
        end select
    end do
    call water%calc_properties_batch(t,p,rho=rho,cp=cp,drho_dT=dt,drho_dP=dp, &
        cv=cv,w=w,u=u,s=s,h=h,nu=nu,p=press,status=status,region_id=regions)
    do i=1,n
        if(status(i)/=0.or.regions(i)/=water%get_region(t(i),p(i))) error stop 'batch dispatch'
        call water%calc_properties(t(i),p(i),ref)
        call close(rho(i),ref%rho);call close(cp(i),ref%cp)
        call close(dt(i),ref%drho_dT);call close(dp(i),ref%drho_dP)
        call close(cv(i),ref%cv);call close(w(i),ref%w)
        call close(u(i),ref%u);call close(s(i),ref%s);call close(h(i),ref%h)
        call close(nu(i),1d0/ref%rho);call close(press(i),p(i))
    end do
    do k=1,3
        do i=1,n
            select case(k)
            case(1)
                t(i)=280d0+mod(i*97,131)*2.5d0;p(i)=30d6+mod(i*67,131)*5d5
            case(2)
                t(i)=700d0+mod(i*97,131)*2d0;p(i)=1d5+mod(i*67,131)*5d4
            case(3)
                t(i)=1200d0+mod(i*97,131)*7d0;p(i)=10d6+mod(i*67,131)*2d5
            end select
        end do
        do m=1,size(sizes)
            call water%calc_properties_batch(t(:sizes(m)),p(:sizes(m)), &
                rho=rho(:sizes(m)),cp=cp(:sizes(m)),drho_dT=dt(:sizes(m)),drho_dP=dp(:sizes(m)))
            do i=1,sizes(m)
                call water%calc_properties(t(i),p(i),ref)
                call close(rho(i),ref%rho);call close(cp(i),ref%cp)
                call close(dt(i),ref%drho_dT);call close(dp(i),ref%drho_dP)
            end do
        end do
    end do
    ! Refresh reference arrays after the homogeneous SIMD stress cases.
    call water%calc_properties_batch(t,p,rho=rho,cp=cp)
    ! Noncontiguous output, partial outputs, pure wrapper, zero-length arrays.
    scratch=-1d0
    call water%calc_properties_batch(t,p,cp=scratch(1:2*n:2))
    if(any(scratch(2:2*n:2)/=-1d0)) error stop 'batch output stride'
    do i=1,n
        call close(scratch(2*i-1),cp(i))
    end do
    t_store=ieee_value(0d0,ieee_quiet_nan);p_store=t_store
    t_store(1:2*n:2)=t;p_store(1:2*n:2)=p
    call water%calc_properties_batch(t_store(1:2*n:2),p_store(1:2*n:2),cp=rho)
    do i=1,n
        call close(rho(i),cp(i))
    end do
    meta_store=-99;region_store=-99
    call water%calc_properties_batch(t,p,status=meta_store(1:2*n:2),region_id=region_store(1:2*n:2))
    if(any(meta_store(1:2*n:2)/=0).or.any(meta_store(2:2*n:2)/=-99)) error stop 'status-only stride'
    do i=1,n
        if(region_store(2*i-1)/=water%get_region(t(i),p(i))) error stop 'region-only dispatch'
    end do
    call evaluate_pure(water,t,p,rho)
    call water%calc_properties_batch(t(:0),p(:0),rho=rho(:0),status=status(:0))
    ! Bad elements must not abort or contaminate neighboring valid inputs.
    nan=ieee_value(0d0,ieee_quiet_nan)
    t(2)=nan;p(4)=-1d0;t(7)=3000d0
    call water%calc_properties_batch(t,p,rho=rho,cp=cp,status=status,region_id=regions)
    do i=1,n
        if(any(i==[2,4,7])) then
            if(status(i)/=-1.or.regions(i)/=-1) error stop 'invalid batch status'
            if(.not.ieee_is_nan(rho(i)).or..not.ieee_is_nan(cp(i))) error stop 'invalid batch output'
        else
            if(status(i)/=0.or..not.ieee_is_finite(cp(i))) error stop 'valid batch neighbor'
        end if
    end do
    ! Native (T,rho) and (T,p) interfaces, every exposed quantity and lane.
    do k=1,2
        do i=1,n
            if(k==1) then
                t(i)=280d0+mod(i,60)*2d0;p(i)=980d0+mod(i,30)
            else
                t(i)=220d0+mod(i,50);p(i)=1d5+mod(i,100)*1d5
            end if
        end do
        if(k==1) then
            call fluid%calc_properties_batch(t,p,rho=rho,cp=cp,drho_dT=dt,drho_dP=dp, &
                cv=cv,w=w,u=u,s=s,h=h,nu=nu,p=press,status=status)
        else
            call ice%calc_properties_batch(t,p,rho=rho,cp=cp,drho_dT=dt,drho_dP=dp, &
                cv=cv,w=w,u=u,s=s,h=h,nu=nu,p=press,status=status)
        end if
        do i=1,n
            if(status(i)/=0) error stop 'valid native batch status'
            if(k==1) then
                call fluid%calc_properties(t(i),p(i),ref)
                call close(press(i),ref%p)
            else
                call ice%calc_properties(t(i),p(i),ref)
                call close(press(i),p(i))
            end if
            call close(rho(i),ref%rho);call close(cp(i),ref%cp)
            call close(dt(i),ref%drho_dT);call close(dp(i),ref%drho_dP)
            call close(cv(i),ref%cv);call close(w(i),ref%w)
            call close(u(i),ref%u);call close(s(i),ref%s);call close(h(i),ref%h)
            call close(nu(i),1d0/ref%rho)
        end do
    end do
    inf=ieee_value(0d0,ieee_positive_inf)
    t=300d0;p=996d0
    t(:3)=[0d0,-1d0,nan];p(4:6)=[0d0,-1d0,inf]
    call fluid%calc_properties_batch(t,p,rho=rho,cp=cp,status=status)
    if(any(status(:6)/=-1).or.any(.not.ieee_is_nan(rho(:6)))) error stop 'invalid native density/temperature'
    if(any(status(7:)/=0).or.any(.not.ieee_is_finite(cp(7:)))) error stop 'valid native neighbors'
    t=260d0;p=1d5;t(3)=nan;p(4)=inf
    call ice%calc_properties_batch(t,p,cp=cp,status=status)
    if(any(status(3:4)/=-1).or.any(.not.ieee_is_nan(cp(3:4)))) error stop 'invalid ice batch'
    ! Ice retains its documented finite extrapolation policy.
    t(:2)=[-1d0,0d0];p(:2)=[-100d0,0d0]
    call ice%calc_properties_batch(t(:2),p(:2),cp=cp(:2),status=status(:2))
    do i=1,2
        if(status(i)/=0) error stop 'ice finite extrapolation policy'
        call ice%calc_cp(t(i),p(i),ref%cp)
        call close(cp(i),ref%cp)
    end do
#ifdef _OPENMP
    call check_parallel_batch()
#endif
contains
    pure subroutine evaluate_pure(model,t,p,rho)
        type(type_iapws97),intent(in)::model
        real(real64),intent(in)::t(:),p(:)
        real(real64),intent(out)::rho(:)
        call model%calc_properties_batch(t,p,rho=rho)
    end subroutine
    subroutine close(actual,expected)
        real(real64),intent(in)::actual,expected
        if(.not.ieee_is_finite(actual).or..not.ieee_is_finite(expected)) error stop 'nonfinite batch result'
        if(abs(actual-expected)>5d-12*max(abs(actual),abs(expected),1d-10)) then
            print *,actual,expected
            error stop 'batch scalar mismatch'
        end if
    end subroutine
end program

#ifdef _OPENMP
subroutine check_parallel_batch
    use iso_fortran_env,only:real64,int64
    use omp_lib,only:omp_get_num_threads
    use iapws,only:type_iapws97
    implicit none
    integer,parameter::n=1027,chunk=257
    type(type_iapws97)::water
    real(real64)::t(n),p(n),serial(4,n),parallel(4,n),ts(n),sat_serial(4,n),sat_parallel(4,n)
    integer::i,first,last,workers
    call water%initialize()
    do i=1,n
        select case(mod(i,4))
        case(0)
            t(i)=280d0+mod(i,100)*3d0;p(i)=30d6
        case(1)
            t(i)=700d0+mod(i,100)*3d0;p(i)=1d5
        case(2)
            t(i)=650d0+mod(i,100)*0.4d0;p(i)=30d6
        case(3)
            t(i)=1200d0+mod(i,100)*8d0;p(i)=10d6
        end select
    end do
    ts=273.15d0+real(mod([(i,i=1,n)],256),real64)*1.45d0
    do first=1,n,chunk
        last=min(n,first+chunk-1)
        call water%calc_properties_batch(t(first:last),p(first:last),rho=serial(1,first:last), &
            cp=serial(2,first:last),drho_dT=serial(3,first:last),drho_dP=serial(4,first:last))
        call water%calc_saturation_properties_batch(ts(first:last),rho_vap=sat_serial(1,first:last), &
            cp_vap=sat_serial(2,first:last),latent_heat=sat_serial(3,first:last),dlatent_heat_dT=sat_serial(4,first:last))
    end do
    workers=0
    !$omp parallel default(none) shared(water,t,p,parallel,workers,ts,sat_parallel) private(first,last) num_threads(4)
    !$omp single
    workers=omp_get_num_threads()
    !$omp end single
    !$omp do schedule(dynamic,1)
    do first=1,n,chunk
        last=min(n,first+chunk-1)
        call water%calc_properties_batch(t(first:last),p(first:last),rho=parallel(1,first:last), &
            cp=parallel(2,first:last),drho_dT=parallel(3,first:last),drho_dP=parallel(4,first:last))
        call water%calc_saturation_properties_batch(ts(first:last),rho_vap=sat_parallel(1,first:last), &
            cp_vap=sat_parallel(2,first:last),latent_heat=sat_parallel(3,first:last), &
            dlatent_heat_dT=sat_parallel(4,first:last))
    end do
    !$omp end do
    !$omp end parallel
    if(workers<2) error stop 'batch parallel check requires multiple workers'
    if(any(transfer(serial,[0_int64])/=transfer(parallel,[0_int64]))) error stop 'parallel batch mismatch'
    if(any(transfer(sat_serial,[0_int64])/=transfer(sat_parallel,[0_int64]))) error stop 'parallel saturation mismatch'
    print *, 'Region queues and SIMD batches agree across workers:',workers
end subroutine check_parallel_batch
#endif
