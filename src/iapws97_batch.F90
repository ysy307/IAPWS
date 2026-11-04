submodule(module_iapws97) iapws97_batch
    use module_iapws_batch
    implicit none
contains
    module pure recursive subroutine calc_properties_batch_iapws97(self,T_in,P_in,rho, drho_dT, drho_dP, cp, cv, w, u, &
            s, h, nu, p, status, region_id)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(out),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer,intent(out),optional::status(:),region_id(:)
        integer::indices(IAPWS_BATCH_TILE,5),counts(5),j,r
        logical::want_values,scalar_model
        call batch_check_shapes(size(T_in),size(P_in),rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p, &
            status, region_id)
        want_values=any([present(rho), present(drho_dT), present(drho_dP), present(cp), present(cv), &
            present(w), present(u), present(s), present(h), present(nu), present(p)])
        select type(self)
        type is(type_iapws97)
            scalar_model=.false.
        class default
            scalar_model=.true.
        end select
        counts=0
        ! Bounded region queues fill SIMD tiles even for interleaved input.
        ! No n-sized permutation, shared workspace, or internal threads.
        do j=1,size(T_in)
            r=self%get_region(T_in(j),P_in(j))
            if(present(region_id)) region_id(j)=r
            if(r<1.or.r>5.or.r==4) then
                if(present(status)) status(j)=IAPWS_BATCH_INVALID
                call batch_invalid(j,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
                cycle
            end if
            if(present(status)) status(j)=IAPWS_BATCH_OK
            if(.not.want_values) cycle
            if(scalar_model) then
                ! A subtype may override calc_properties. Preserve that
                ! behavior rather than bypassing it through private kernels.
                call evaluate_region(self,3,[j],T_in,P_in,rho,drho_dT,drho_dP,cp,cv,w,u,s,h,nu,p)
                cycle
            end if
            counts(r)=counts(r)+1;indices(counts(r),r)=j
            if(counts(r)==IAPWS_BATCH_TILE) then
                call evaluate_region(self,r,indices(:,r),T_in,P_in,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
                counts(r)=0
            end if
        end do
        do r=1,5
            if(counts(r)==0) cycle
            call evaluate_region(self,r,indices(:counts(r),r),T_in,P_in,rho, drho_dT, drho_dP, cp, cv, w, u, &
            s, h, nu, p)
        end do
    end subroutine

    pure recursive subroutine evaluate_region(self,region,indices,T_in,P_in,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
        class(type_iapws97),intent(in)::self
        integer,intent(in)::region,indices(:)
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(inout),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        real(real64)::t(IAPWS_BATCH_TILE),pressure(IAPWS_BATCH_TILE),g(IAPWS_BATCH_TILE,6),values(1,11)
        type(type_iapws_property)::props
        integer::k,id,count
        count=size(indices)
        do k=1,count
            id=indices(k);t(k)=T_in(id);pressure(k)=P_in(id)
        end do
        select case(region)
        case(1)
            call self%region1%calc_gamma_batch(t(:count),pressure(:count),g(:count,:))
        case(2)
            call self%region2%calc_gamma_batch(t(:count),pressure(:count),g(:count,:))
        case(5)
            call self%region5%calc_gamma_batch(t(:count),pressure(:count),g(:count,:))
        case(3)
            ! Preserve each lane's exact root-solving history and safeguards.
            do k=1,count
                id=indices(k)
                call self%calc_properties(t(k),pressure(k),props)
                values(1,:)=[props%rho,props%drho_dT,props%drho_dP,props%cp,props%cv, &
                    props%w,props%u,props%s,props%h,1d0/props%rho,pressure(k)]
                call batch_values_outputs([id],values,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
            end do
            return
        end select
        call batch_gibbs_outputs(indices,t(:count),pressure(:count),g(:count,:),rho, drho_dT, drho_dP, cp, cv, &
            w, u, s, h, nu, p)
    end subroutine
end submodule iapws97_batch
