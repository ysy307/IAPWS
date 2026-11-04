submodule(module_iapws97) iapws97_base
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    !> Invalid region identifier
    integer(int32), parameter :: IAPWS97_INVALID = -1
    !> Region 1 (Liquid water) identifier
    integer(int32), parameter :: IAPWS97_REGION_1 = 1
    !> Region 2 (Steam) identifier
    integer(int32), parameter :: IAPWS97_REGION_2 = 2
    !> Region 3 (Critical / Supercritical) identifier
    integer(int32), parameter :: IAPWS97_REGION_3 = 3
    !> Region 5 (High temperature) identifier
    integer(int32), parameter :: IAPWS97_REGION_5 = 5

    !> Global maximum temperature limit [K]
    real(real64), parameter :: IAPWS97_LIMIT_T_MAX = 2273.15d0
    !> Global maximum pressure limit [Pa]
    real(real64), parameter :: IAPWS97_LIMIT_P_MAX = 100.0d6

    !> Region 5 minimum temperature [K]
    real(real64), parameter :: IAPWS97_R5_T_MIN = 1073.15d0
    !> Region 5 maximum pressure [Pa]
    real(real64), parameter :: IAPWS97_R5_P_MAX = 50.0d6
    ! Phase data shared by saturation values, partials, and total derivatives.
    type :: saturation_phase
        real(real64) :: rho, cp, h, rho_T, rho_p, cp_T, cp_p, h_p
    end type
contains

    !> Initialize the IAPWS-97 model instance.
    !> Initializes all sub-region handlers and auxiliary equations.
    module pure elemental subroutine initialize_type_iapws97(self)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(inout) :: self

        call self%auxiliary%initialize()
        call self%region1%initialize()
        call self%region2%initialize()
        call self%region3%initialize()
        call self%region4%initialize()
        call self%region5%initialize()

    end subroutine initialize_type_iapws97

    !> Determine the IAPWS-97 region for a given temperature and pressure.
    !> Checks validity limits and boundaries between regions (Saturation, B23).

    module pure elemental function get_region_iapws97(self, T_in, p_in) result(region_id)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Region ID result
        integer(int32) :: region_id

        real(real64) :: p_boundary

        ! 1. Global Validity Check
        ! Region 1, 2, 3: p <= 100 MPa
        ! Region 5: p <= 50 MPa -> Checked individually
        ! T >= 243.15 K
        region_id = IAPWS97_INVALID
        if (.not.ieee_is_finite(T_in).or..not.ieee_is_finite(p_in)) return
        if (p_in <= 0.0d0 .or. p_in > 100.0d6 .or. &
            T_in < 243.15d0 .or. T_in > 2273.15d0) then
            region_id = IAPWS97_INVALID
            return
        end if

        ! 2. Hierarchical determination by temperature (Based on Fig. 1)

        ! --- High Temperature (Region 5) ---
        ! 1073.15 K < T <= 2273.15 K
        if (T_in > 1073.15d0) then
            if (p_in <= 50.0d6) then
                region_id = IAPWS97_REGION_5
            else
                region_id = IAPWS97_INVALID ! Region 5 max P is 50MPa
            end if
            return
        end if

        ! --- Region 2 Only Area ---
        ! 863.15 K < T <= 1073.15 K
        ! This temp range is Region 2 for all pressures up to 100MPa
        if (T_in > 863.15d0) then
            region_id = IAPWS97_REGION_2
            return
        end if

        ! --- Boundary between Region 2 and Region 3 ---
        ! 623.15 K < T <= 863.15 K
        ! Boundary: B23 equation (Eq. 5)
        if (T_in > 623.15d0) then
            p_boundary = self%auxiliary%calc_p_boundary(T_in)

            if (p_in > p_boundary) then
                region_id = IAPWS97_REGION_3
            else
                region_id = IAPWS97_REGION_2
            end if
            return
        end if

        ! --- Boundary between Region 1 and Region 2 (Saturation Curve) ---
        ! 273.15 K <= T <= 623.15 K
        ! Boundary: Saturation pressure equation (Eq. 30)
        p_boundary = self%region4%calc_psat(T_in)

        if (p_in > p_boundary) then
            region_id = IAPWS97_REGION_1 ! Compressed liquid
        else
            region_id = IAPWS97_REGION_2 ! Superheated (or metastable) steam
        end if

    end function get_region_iapws97

    !> Calculate all thermodynamic properties for a given state.
    !> Dispatches calculation to the appropriate region based on \( T \) and \( P \).
    module pure elemental subroutine calc_properties_iapws97(self, T_in, p_in, property)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Output property structure
        type(type_iapws_property), intent(inout) :: property

        property%region_id = self%get_region(T_in, p_in)
        select case (property%region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_2)
            call self%region2%calc_properties(T_in, p_in, property)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, property%rho)
            call self%region3%calc_properties(T_in, property%rho, property)
        case (IAPWS97_REGION_5)
            call self%region5%calc_properties(T_in, p_in, property)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_properties_iapws97."
        end select

    end subroutine calc_properties_iapws97

    !> Calculate specific volume \( \nu \).
    module pure elemental subroutine calc_nu_iapws97(self, T_in, p_in, nu)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific volume [m^3/kg]
        real(real64), intent(inout) :: nu

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_nu(T_in, p_in, nu)
        case (IAPWS97_REGION_2)
            call self%region2%calc_nu(T_in, p_in, nu)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            nu = 1.0d0 / rho
        case (IAPWS97_REGION_5)
            call self%region5%calc_nu(T_in, p_in, nu)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_nu_iapws97."
        end select
    end subroutine calc_nu_iapws97

    !> Calculate density \( \rho \).
    module pure elemental subroutine calc_rho_iapws97(self, T_in, p_in, rho)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Density [kg/m^3]
        real(real64), intent(inout) :: rho

        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_rho(T_in, p_in, rho)
        case (IAPWS97_REGION_2)
            call self%region2%calc_rho(T_in, p_in, rho)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
        case (IAPWS97_REGION_5)
            call self%region5%calc_rho(T_in, p_in, rho)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_rho_iapws97."
        end select

    end subroutine calc_rho_iapws97

    module pure elemental subroutine calc_drho_dT_iapws97(self, T_in, p_in, drho_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dT

        real(real64) :: rho
        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dT(T_in, p_in, drho_dT)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_drho_dT(T_in, rho, drho_dT)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dT(T_in, p_in, drho_dT)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_drho_dT_iapws97."
        end select

    end subroutine calc_drho_dT_iapws97

    module pure elemental subroutine calc_drho_dp_iapws97(self, T_in, p_in, drho_dp)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in
        real(real64), intent(in) :: p_in
        real(real64), intent(inout) :: drho_dp

        real(real64) :: rho
        integer(int32) :: region_id
        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_2)
            call self%region2%calc_drho_dp(T_in, p_in, drho_dp)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_drho_dp(T_in, rho, drho_dp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_drho_dp(T_in, p_in, drho_dp)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_drho_dp_iapws97."
        end select

    end subroutine calc_drho_dp_iapws97

    !> Calculate specific internal energy \( u \).
    module pure elemental subroutine calc_u_iapws97(self, T_in, p_in, u)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific internal energy [J/kg]
        real(real64), intent(inout) :: u

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_u(T_in, p_in, u)
        case (IAPWS97_REGION_2)
            call self%region2%calc_u(T_in, p_in, u)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_u(T_in, rho, u)
        case (IAPWS97_REGION_5)
            call self%region5%calc_u(T_in, p_in, u)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_u_iapws97."
        end select

    end subroutine calc_u_iapws97

    !> Calculate specific enthalpy \( h \).
    module pure elemental subroutine calc_h_iapws97(self, T_in, p_in, h)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific enthalpy [J/kg]
        real(real64), intent(inout) :: h

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_h(T_in, p_in, h)
        case (IAPWS97_REGION_2)
            call self%region2%calc_h(T_in, p_in, h)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_h(T_in, rho, h)
        case (IAPWS97_REGION_5)
            call self%region5%calc_h(T_in, p_in, h)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_h_iapws97."
        end select

    end subroutine calc_h_iapws97

    !> Calculate specific entropy \( s \).
    module pure elemental subroutine calc_s_iapws97(self, T_in, p_in, s)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific entropy [J/(kg K)]
        real(real64), intent(inout) :: s

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_s(T_in, p_in, s)
        case (IAPWS97_REGION_2)
            call self%region2%calc_s(T_in, p_in, s)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_s(T_in, rho, s)
        case (IAPWS97_REGION_5)
            call self%region5%calc_s(T_in, p_in, s)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_s_iapws97."
        end select

    end subroutine calc_s_iapws97

    !> Calculate specific isobaric heat capacity \( c_p \).
    module pure elemental subroutine calc_cp_iapws97(self, T_in, p_in, cp)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific isobaric heat capacity [J/(kg K)]
        real(real64), intent(inout) :: cp

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)
        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_cp(T_in, p_in, cp)
        case (IAPWS97_REGION_2)
            call self%region2%calc_cp(T_in, p_in, cp)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_cp(T_in, rho, cp)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cp(T_in, p_in, cp)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_cp_iapws97."
        end select

    end subroutine calc_cp_iapws97

    !> Calculate specific isochoric heat capacity \( c_v \).
    module pure elemental subroutine calc_cv_iapws97(self, T_in, p_in, cv)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Specific isochoric heat capacity [J/(kg K)]
        real(real64), intent(inout) :: cv

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)
        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_cv(T_in, p_in, cv)
        case (IAPWS97_REGION_2)
            call self%region2%calc_cv(T_in, p_in, cv)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_cv(T_in, rho, cv)
        case (IAPWS97_REGION_5)
            call self%region5%calc_cv(T_in, p_in, cv)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_cv_iapws97."
        end select

    end subroutine calc_cv_iapws97

    !> Calculate speed of sound \( w \).
    module pure elemental subroutine calc_w_iapws97(self, T_in, p_in, w)
        implicit none
        !> IAPWS-97 instance
        class(type_iapws97), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: T_in
        !> Pressure [Pa]
        real(real64), intent(in) :: p_in
        !> Speed of sound [m/s]
        real(real64), intent(inout) :: w

        integer(int32) :: region_id
        real(real64) :: rho

        region_id = self%get_region(T_in, p_in)

        select case (region_id)
        case (IAPWS97_REGION_1)
            call self%region1%calc_w(T_in, p_in, w)
        case (IAPWS97_REGION_2)
            call self%region2%calc_w(T_in, p_in, w)
        case (IAPWS97_REGION_3)
            call region3_density(self, T_in, p_in, rho)
            call self%region3%calc_w(T_in, rho, w)
        case (IAPWS97_REGION_5)
            call self%region5%calc_w(T_in, p_in, w)
        case default
            error stop "IAPWS97: State out of bounds (T or P invalid) in calc_w_iapws97."
        end select

    end subroutine calc_w_iapws97

    !> cp temperature derivative at fixed pressure [J/(kg K^2)].
    module pure elemental subroutine calc_dcp_dT_iapws97(self, T_in, p_in, dcp_dT)
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dT
        real(real64) :: unused
        call cp_derivatives(self,T_in,p_in,dcp_dT,unused)
    end subroutine

    !> cp pressure derivative at fixed temperature [J/(kg K Pa)].
    module pure elemental subroutine calc_dcp_dp_iapws97(self, T_in, p_in, dcp_dp)
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in, p_in
        real(real64), intent(out) :: dcp_dp
        real(real64) :: unused
        call cp_derivatives(self,T_in,p_in,unused,dcp_dp)
    end subroutine

    pure elemental subroutine cp_derivatives(self,T,p,ct,cp)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T,p
        real(real64),intent(out)::ct,cp
        real(real64)::rho
        ct=ieee_value(0d0,ieee_quiet_nan); cp=ct
        if (.not.ieee_is_finite(T).or..not.ieee_is_finite(p)) return
        if (T==647.096d0.and.abs(p-22.064d6)<1d-3) return
        select case(self%get_region(T,p))
        case(1)
            call self%region1%calc_cp_derivatives(T,p,ct,cp)
        case(2)
            call self%region2%calc_cp_derivatives(T,p,ct,cp)
        case(3)
            call region3_density(self,T,p,rho)
            call self%region3%calc_cp_derivatives(T,rho,ct,cp)
        case(5)
            call self%region5%calc_cp_derivatives(T,p,ct,cp)
        end select
    end subroutine

    !> Safeguarded region-3 solve shared by values and ordinary derivatives.
    pure elemental subroutine region3_density(self,T,p,rho)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T,p
        real(real64),intent(out)::rho
        real(real64)::lo,hi,mid,pm,slope,trial
        integer::i
        if(T<self%region3%T_c) then
            call saturation_root(self,T,p,p<=self%region4%calc_psat(T),rho)
            return
        end if
        lo=1d0; hi=800d0; mid=(lo+hi)/2d0
        do i=1,100
            call self%region3%calc_p(T,mid,pm)
            if(pm>p) then
                hi=mid
            else
                lo=mid
            end if
            if(hi-lo<2d-11) exit
            call self%region3%calc_p_rho(T,mid,slope)
            trial=(lo+hi)/2d0
            if(slope>0d0) then
                if(mid-(pm-p)/slope>lo.and.mid-(pm-p)/slope<hi) trial=mid-(pm-p)/slope
            end if
            mid=trial
        end do
        rho=(lo+hi)/2d0
        call self%region3%calc_p(T,rho,pm)
        if(abs(pm-p)>0.01d0) rho=ieee_value(0d0,ieee_quiet_nan)
    end subroutine

    !> Stable phase root, bracketed by spinodals to prevent phase switching.
    pure elemental subroutine saturation_root(self,T,p,vapor,rho)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T,p
        logical,intent(in)::vapor
        real(real64),intent(out)::rho
        real(real64)::lo,hi,mid,slope,plo,phi,pm,trial
        integer::i
        rho=ieee_value(0d0,ieee_quiet_nan)
        ! Locate the spinodal on the requested side of the critical density.
        if(vapor) then
            lo=1d0; hi=self%region3%rho_c
        else
            lo=self%region3%rho_c; hi=800d0
        end if
        do i=1,70
            mid=(lo+hi)/2d0
            call self%region3%calc_p_rho(T,mid,slope)
            if ((slope>0d0).eqv.vapor) then
                lo=mid
            else
                hi=mid
            end if
        end do
        if(vapor) then
            hi=(lo+hi)/2d0; lo=1d0
        else
            lo=(lo+hi)/2d0; hi=800d0
        end if
        call self%region3%calc_p(T,lo,plo)
        call self%region3%calc_p(T,hi,phi)
        if(p<plo.or.p>phi) return
        mid=(lo+hi)/2d0
        do i=1,100
            call self%region3%calc_p(T,mid,pm)
            if(pm>p) then
                hi=mid
            else
                lo=mid
            end if
            if(hi-lo<2d-11) exit
            call self%region3%calc_p_rho(T,mid,slope)
            trial=(lo+hi)/2d0
            if(slope>0d0) then
                if(mid-(pm-p)/slope>lo.and.mid-(pm-p)/slope<hi) trial=mid-(pm-p)/slope
            end if
            mid=trial
        end do
        rho=(lo+hi)/2d0
        call self%region3%calc_p(T,rho,pm)
        call self%region3%calc_p_rho(T,rho,slope)
        if(abs(pm-p)>0.01d0.or.slope<=0d0) rho=ieee_value(0d0,ieee_quiet_nan)
    end subroutine

    pure elemental subroutine saturation_state(self,T,vap,liq,slope,derivatives)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T
        type(saturation_phase),intent(out)::vap,liq
        real(real64),intent(out)::slope
        logical,intent(in)::derivatives
        real(real64)::p,nan
        nan=ieee_value(0d0,ieee_quiet_nan)
        vap=saturation_phase(nan,nan,nan,nan,nan,nan,nan,nan); liq=vap; slope=nan
        if(.not.ieee_is_finite(T)) return
        if(T<273.15d0.or.T>647.096d0) return
        if(T==647.096d0) then
            if(.not.derivatives) then
                vap%rho=322d0; liq%rho=322d0
                call self%region3%calc_h(T,322d0,vap%h)
                liq%h=vap%h
            end if
            return
        end if
        p=self%region4%calc_psat(T)
        slope=self%region4%calc_dpsat_dT(T)
        if(T<=623.15d0) then
            call saturation_gibbs(self%region2,T,p,vap)
            call saturation_gibbs(self%region1,T,p,liq)
            if(derivatives) then
                call self%region2%calc_cp_derivatives(T,p,vap%cp_T,vap%cp_p)
                call self%region1%calc_cp_derivatives(T,p,liq%cp_T,liq%cp_p)
            end if
        else
            call saturation_root(self,T,p,.true.,vap%rho)
            call saturation_root(self,T,p,.false.,liq%rho)
            if(.not.ieee_is_finite(vap%rho).or..not.ieee_is_finite(liq%rho)) then
                ! The independent IF97 fits can lose a stable root extremely
                ! close to Tc. Reject the whole pair instead of leaking one phase.
                vap%rho=nan
                liq%rho=nan
                return
            end if
            call saturation_helmholtz(self%region3,T,vap,derivatives)
            call saturation_helmholtz(self%region3,T,liq,derivatives)
        end if
    end subroutine

    pure elemental subroutine saturation_gibbs(region,T,p,phase)
        class(abst_iapws_gibbs),intent(in)::region
        real(real64),intent(in)::T,p
        type(saturation_phase),intent(inout)::phase
        type(type_iapws_gibbs_coefficient)::g
        call region%calc_gamma(T,p,g)
        phase%rho=1d0/g%g_p
        phase%cp=-T*g%g_tt
        phase%h=g%g-T*g%g_t
        phase%rho_T=-g%g_pt/g%g_p**2
        phase%rho_p=-g%g_pp/g%g_p**2
        phase%h_p=g%g_p-T*g%g_pt
    end subroutine

    pure elemental subroutine saturation_helmholtz(region,T,phase,derivatives)
        class(type_iapws97_region3),intent(in)::region
        real(real64),intent(in)::T
        type(saturation_phase),intent(inout)::phase
        logical,intent(in)::derivatives
        type(type_iapws_helmholtz_property)::f
        real(real64)::pt,pr
        call region%calc_phi(region%T_c/T,phase%rho/region%rho_c,f)
        call region%calc_cp(T,phase%rho,phase%cp,f)
        call region%calc_h(T,phase%rho,phase%h,f)
        call region%calc_p_T(T,phase%rho,pt,f)
        call region%calc_p_rho(T,phase%rho,pr,f)
        phase%rho_T=-pt/pr
        phase%rho_p=1d0/pr
        ! h_p|T = v - T*v_T|p, valid also for the Helmholtz EOS.
        phase%h_p=1d0/phase%rho+T*phase%rho_T/phase%rho**2
        if(derivatives) call region%calc_cp_derivatives(T,phase%rho,phase%cp_T,phase%cp_p)
    end subroutine

    module pure elemental subroutine calc_latent_heat_iapws97(self, latent_heat, T_in, p_in)
        class(type_iapws97), intent(in) :: self
        real(real64), intent(inout) :: latent_heat
        real(real64), intent(in), optional :: T_in, p_in
        type(saturation_phase)::vap,liq
        real(real64)::T,slope
        latent_heat=ieee_value(0d0,ieee_quiet_nan)
        if(present(T_in).eqv.present(p_in)) return
        if(present(T_in)) then
            T=T_in
        else
            if(.not.ieee_is_finite(p_in)) return
            if(p_in<self%region4%calc_psat(273.15d0).or.p_in>22.064d6) return
            T=self%region4%calc_tsat(p_in)
            if(p_in==22.064d6) T=647.096d0
        end if
        call saturation_state(self,T,vap,liq,slope,.false.)
        latent_heat=vap%h-liq%h
    end subroutine

    module pure elemental subroutine calc_latent_heat_dT_iapws97(self,T_in,dlatent_heat_dT)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(out)::dlatent_heat_dT
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.true.)
        dlatent_heat_dT=(vap%cp-liq%cp)+(vap%h_p-liq%h_p)*slope
    end subroutine

    module pure elemental subroutine calc_saturation_density_iapws97(self,T_in,rho_vap,rho_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(inout),optional::rho_vap,rho_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.false.)
        if(present(rho_vap)) rho_vap=vap%rho
        if(present(rho_liq)) rho_liq=liq%rho
    end subroutine

    module pure elemental subroutine calc_saturation_cp_iapws97(self,T_in,cp_vap,cp_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(inout),optional::cp_vap,cp_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.false.)
        if(present(cp_vap)) cp_vap=vap%cp
        if(present(cp_liq)) cp_liq=liq%cp
    end subroutine

    module pure elemental subroutine calc_saturation_drho_dT_iapws97(self,T_in,drho_dT_vap,drho_dT_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(inout),optional::drho_dT_vap,drho_dT_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.true.)
        if(present(drho_dT_vap)) drho_dT_vap=vap%rho_T
        if(present(drho_dT_liq)) drho_dT_liq=liq%rho_T
    end subroutine

    module pure elemental subroutine calc_saturation_drho_dP_iapws97(self,T_in,drho_dP_vap,drho_dP_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(inout),optional::drho_dP_vap,drho_dP_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.true.)
        if(present(drho_dP_vap)) drho_dP_vap=vap%rho_p
        if(present(drho_dP_liq)) drho_dP_liq=liq%rho_p
    end subroutine

    module pure elemental subroutine calc_saturation_dcp_dT_iapws97(self,T_in,dcp_dT_vap,dcp_dT_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(out),optional::dcp_dT_vap,dcp_dT_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.true.)
        if(present(dcp_dT_vap)) dcp_dT_vap=vap%cp_T+vap%cp_p*slope
        if(present(dcp_dT_liq)) dcp_dT_liq=liq%cp_T+liq%cp_p*slope
    end subroutine

    module pure elemental subroutine calc_saturation_density_dT_iapws97(self,T_in,drho_dT_vap,drho_dT_liq)
        class(type_iapws97),intent(in)::self
        real(real64),intent(in)::T_in
        real(real64),intent(out),optional::drho_dT_vap,drho_dT_liq
        type(saturation_phase)::vap,liq
        real(real64)::slope
        call saturation_state(self,T_in,vap,liq,slope,.true.)
        if(present(drho_dT_vap)) drho_dT_vap=vap%rho_T+vap%rho_p*slope
        if(present(drho_dT_liq)) drho_dT_liq=liq%rho_T+liq%rho_p*slope
    end subroutine

end submodule iapws97_base
