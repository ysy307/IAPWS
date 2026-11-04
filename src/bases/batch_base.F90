submodule(module_iapws) iapws_batch_base
    use ieee_arithmetic,only:ieee_is_finite
    use module_iapws_batch
    implicit none
contains
    module pure recursive subroutine calc_gamma_batch_gibbs(self,T_in,P_in,coef)
        class(abst_iapws_gibbs),intent(in)::self
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(out)::coef(:,:)
        type(type_iapws_gibbs_coefficient)::g
        integer::j
        if(size(P_in)/=size(T_in).or.size(coef,1)/=size(T_in).or.size(coef,2)/=6) &
            error stop 'IAPWS batch: Gibbs coefficient shape'
        do j=1,size(T_in)
            call self%calc_gamma(T_in(j),P_in(j),g)
            coef(j,:)=[g%g,g%g_p,g%g_t,g%g_pp,g%g_tt,g%g_pt]
        end do
    end subroutine

    module pure recursive subroutine calc_properties_batch_gibbs(self,T_in,P_in,rho, drho_dT, drho_dP, cp, cv, &
            w, u, s, h, nu, p, status, region_id)
        class(abst_iapws_gibbs),intent(in)::self
        real(real64),intent(in)::T_in(:),P_in(:)
        real(real64),intent(out),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer,intent(out),optional::status(:),region_id(:)
        type(type_iapws_gibbs_coefficient)::g
        integer::j
        logical::want_gamma
        call batch_check_shapes(size(T_in),size(P_in),rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p, &
            status, region_id)
        want_gamma=any([present(rho), present(drho_dT), present(drho_dP), present(cp), present(cv), &
            present(w), present(u), present(s), present(h), present(nu)])
        ! Conservative inherited path: preserve overridden property bindings
        ! and the full coefficient object of arbitrary external model classes.
        ! Concrete IF97 kernels use the direct SIMD projection separately.
        do j=1,size(T_in)
            if(.not.ieee_is_finite(T_in(j)).or..not.ieee_is_finite(P_in(j))) then
                if(present(status)) status(j)=IAPWS_BATCH_INVALID
                if(present(region_id)) region_id(j)=-1
                call batch_invalid(j,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
                cycle
            end if
            if(present(status)) status(j)=IAPWS_BATCH_OK
            if(present(region_id)) region_id(j)=0
            if(present(p)) p(j)=P_in(j)
            if(.not.want_gamma) cycle
            call self%calc_gamma(T_in(j),P_in(j),g)
            if(present(rho)) call self%calc_rho(T_in(j),P_in(j),rho(j),g)
            if(present(drho_dT)) call self%calc_drho_dT(T_in(j),P_in(j),drho_dT(j),g)
            if(present(drho_dP)) call self%calc_drho_dP(T_in(j),P_in(j),drho_dP(j),g)
            if(present(cp)) call self%calc_cp(T_in(j),P_in(j),cp(j),g)
            if(present(cv)) call self%calc_cv(T_in(j),P_in(j),cv(j),g)
            if(present(w)) call self%calc_w(T_in(j),P_in(j),w(j),g)
            if(present(u)) call self%calc_u(T_in(j),P_in(j),u(j),g)
            if(present(s)) call self%calc_s(T_in(j),P_in(j),s(j),g)
            if(present(h)) call self%calc_h(T_in(j),P_in(j),h(j),g)
            if(present(nu)) call self%calc_nu(T_in(j),P_in(j),nu(j),g)
        end do
    end subroutine

    module pure recursive subroutine calc_properties_batch_helmholtz(self,T_in,rho_in,rho, drho_dT, drho_dP, &
            cp, cv, w, &
            u, s, h, nu, p, status, region_id)
        class(abst_iapws_helmholtz),intent(in)::self
        real(real64),intent(in)::T_in(:),rho_in(:)
        real(real64),intent(out),optional::rho(:), drho_dT(:), drho_dP(:), cp(:), cv(:), w(:), u(:), s(:), &
            h(:), nu(:), p(:)
        integer,intent(out),optional::status(:),region_id(:)
        type(type_iapws_helmholtz_property)::f
        integer::j
        logical::valid,want_phi
        call batch_check_shapes(size(T_in),size(rho_in),rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p, &
            status, region_id)
        want_phi=any([present(drho_dT), present(drho_dP), present(cp), present(cv), present(w), present(u), &
            present(s), present(h), present(p)])
        do j=1,size(T_in)
            valid=ieee_is_finite(T_in(j)).and.ieee_is_finite(rho_in(j))
            if(valid) valid=T_in(j)>0d0.and.rho_in(j)>0d0
            if(.not.valid) then
                if(present(status)) status(j)=IAPWS_BATCH_INVALID
                if(present(region_id)) region_id(j)=-1
                call batch_invalid(j,rho, drho_dT, drho_dP, cp, cv, w, u, s, h, nu, p)
                cycle
            end if
            if(present(status)) status(j)=IAPWS_BATCH_OK
            if(present(region_id)) region_id(j)=0
            if(present(rho)) rho(j)=rho_in(j)
            if(present(nu)) nu(j)=1d0/rho_in(j)
            if(.not.want_phi) cycle
            call self%calc_phi(self%T_c/T_in(j),rho_in(j)/self%rho_c,f)
            if(present(drho_dT)) call self%calc_drho_dT(T_in(j),rho_in(j),drho_dT(j),f)
            if(present(drho_dP)) call self%calc_drho_dP(T_in(j),rho_in(j),drho_dP(j),f)
            if(present(cp)) call self%calc_cp(T_in(j),rho_in(j),cp(j),f)
            if(present(cv)) call self%calc_cv(T_in(j),rho_in(j),cv(j),f)
            if(present(w)) call self%calc_w(T_in(j),rho_in(j),w(j),f)
            if(present(u)) call self%calc_u(T_in(j),rho_in(j),u(j),f)
            if(present(s)) call self%calc_s(T_in(j),rho_in(j),s(j),f)
            if(present(h)) call self%calc_h(T_in(j),rho_in(j),h(j),f)
            if(present(p)) call self%calc_p(T_in(j),rho_in(j),p(j),f)
        end do
    end subroutine
end submodule iapws_batch_base
