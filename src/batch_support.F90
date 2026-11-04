module module_iapws_batch
    use iso_fortran_env,only:real64
    use ieee_arithmetic,only:ieee_value,ieee_quiet_nan
    implicit none
    private
#ifndef IAPWS_BATCH_TILE_SIZE
#define IAPWS_BATCH_TILE_SIZE 32
#endif
#if IAPWS_BATCH_TILE_SIZE < 4 || IAPWS_BATCH_TILE_SIZE > 128
#error IAPWS_BATCH_TILE_SIZE_must_be_between_4_and_128
#endif
    integer,parameter,public::IAPWS_BATCH_TILE=IAPWS_BATCH_TILE_SIZE
    integer,parameter,public::IAPWS_BATCH_OK=0,IAPWS_BATCH_INVALID=-1
    public::batch_check_shapes,batch_invalid,batch_gibbs_outputs,batch_values_outputs
contains
    pure recursive subroutine batch_check_shapes(n,second,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p, status, region_id)
        integer,intent(in)::n,second
        real(real64),intent(in),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer,intent(in),optional::status(:),region_id(:)
        if(second/=n) error stop "IAPWS batch: input sizes differ"
        if(present(rho)) then
            if(size(rho)/=n) error stop "IAPWS batch: rho size differs from input"
        end if
        if(present(drho_dT)) then
            if(size(drho_dT)/=n) error stop "IAPWS batch: drho_dT size differs from input"
        end if
        if(present(drho_dP)) then
            if(size(drho_dP)/=n) error stop "IAPWS batch: drho_dP size differs from input"
        end if
        if(present(cp)) then
            if(size(cp)/=n) error stop "IAPWS batch: cp size differs from input"
        end if
        if(present(cv)) then
            if(size(cv)/=n) error stop "IAPWS batch: cv size differs from input"
        end if
        if(present(w)) then
            if(size(w)/=n) error stop "IAPWS batch: w size differs from input"
        end if
        if(present(u)) then
            if(size(u)/=n) error stop "IAPWS batch: u size differs from input"
        end if
        if(present(s)) then
            if(size(s)/=n) error stop "IAPWS batch: s size differs from input"
        end if
        if(present(h)) then
            if(size(h)/=n) error stop "IAPWS batch: h size differs from input"
        end if
        if(present(nu)) then
            if(size(nu)/=n) error stop "IAPWS batch: nu size differs from input"
        end if
        if(present(p)) then
            if(size(p)/=n) error stop "IAPWS batch: p size differs from input"
        end if
        if(present(status)) then
            if(size(status)/=n) error stop "IAPWS batch: status size differs from input"
        end if
        if(present(region_id)) then
            if(size(region_id)/=n) error stop "IAPWS batch: region_id size differs from input"
        end if
    end subroutine

    pure recursive subroutine batch_invalid(index,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
        integer,intent(in)::index
        real(real64),intent(inout),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        real(real64)::nan
        nan=ieee_value(0d0,ieee_quiet_nan)
        if(present(rho)) rho(index)=nan
        if(present(drho_dT)) drho_dT(index)=nan
        if(present(drho_dP)) drho_dP(index)=nan
        if(present(cp)) cp(index)=nan
        if(present(cv)) cv(index)=nan
        if(present(w)) w(index)=nan
        if(present(u)) u(index)=nan
        if(present(s)) s(index)=nan
        if(present(h)) h(index)=nan
        if(present(nu)) nu(index)=nan
        if(present(p)) p(index)=nan
    end subroutine

    pure recursive subroutine batch_gibbs_outputs(indices,t,pressure,g,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
        integer,intent(in)::indices(:)
        real(real64),intent(in)::t(:),pressure(:),g(:,:)
        real(real64),intent(inout),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer::j
        real(real64)::w_sq
        if(present(rho)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                rho(indices(j))=1d0/g(j,2)
            end do
        end if
        if(present(drho_dT)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                drho_dT(indices(j))=-g(j,6)/(g(j,2)**2)
            end do
        end if
        if(present(drho_dP)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                drho_dP(indices(j))=-g(j,4)/(g(j,2)**2)
            end do
        end if
        if(present(cp)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                cp(indices(j))=-t(j)*g(j,5)
            end do
        end if
        if(present(cv)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                cv(indices(j))=-t(j)*g(j,5)+t(j)*(g(j,6)**2)/g(j,4)
            end do
        end if
        if(present(w)) then
#if !defined(__NVCOMPILER)
            !$omp simd private(w_sq)
#endif
            do j=1,size(indices)
                w_sq=(g(j,2)**2)/((g(j,6)**2/g(j,5))-g(j,4))
                w(indices(j))=sqrt(max(w_sq,0d0))
            end do
        end if
        if(present(u)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                u(indices(j))=g(j,1)-t(j)*g(j,3)-pressure(j)*g(j,2)
            end do
        end if
        if(present(s)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                s(indices(j))=-g(j,3)
            end do
        end if
        if(present(h)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                h(indices(j))=g(j,1)-t(j)*g(j,3)
            end do
        end if
        if(present(nu)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                nu(indices(j))=g(j,2)
            end do
        end if
        if(present(p)) then
#if !defined(__NVCOMPILER)
            !$omp simd
#endif
            do j=1,size(indices)
                p(indices(j))=pressure(j)
            end do
        end if
    end subroutine

    pure recursive subroutine batch_values_outputs(indices,values,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
        integer,intent(in)::indices(:)
        real(real64),intent(in)::values(:,:)
        real(real64),intent(inout),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer::j
        if(present(rho)) then
            do j=1,size(indices)
                rho(indices(j))=values(j,1)
            end do
        end if
        if(present(drho_dT)) then
            do j=1,size(indices)
                drho_dT(indices(j))=values(j,2)
            end do
        end if
        if(present(drho_dP)) then
            do j=1,size(indices)
                drho_dP(indices(j))=values(j,3)
            end do
        end if
        if(present(cp)) then
            do j=1,size(indices)
                cp(indices(j))=values(j,4)
            end do
        end if
        if(present(cv)) then
            do j=1,size(indices)
                cv(indices(j))=values(j,5)
            end do
        end if
        if(present(w)) then
            do j=1,size(indices)
                w(indices(j))=values(j,6)
            end do
        end if
        if(present(u)) then
            do j=1,size(indices)
                u(indices(j))=values(j,7)
            end do
        end if
        if(present(s)) then
            do j=1,size(indices)
                s(indices(j))=values(j,8)
            end do
        end if
        if(present(h)) then
            do j=1,size(indices)
                h(indices(j))=values(j,9)
            end do
        end if
        if(present(nu)) then
            do j=1,size(indices)
                nu(indices(j))=values(j,10)
            end do
        end if
        if(present(p)) then
            do j=1,size(indices)
                p(indices(j))=values(j,11)
            end do
        end if
    end subroutine
end module module_iapws_batch
