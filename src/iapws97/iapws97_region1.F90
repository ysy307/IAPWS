submodule(module_iapws97) iapws97_region1
    !------------------------------------------------------------------------------------------
    ! Region1: Saturated liquid water (IAPWS-IF97)
    !------------------------------------------------------------------------------------------
    integer(int32), parameter :: N1_terms = 34
    integer(int32), parameter :: I_r1(N1_terms) = [ &
                                 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
                                 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, &
                                 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, &
                                 29, 30, 31, 32]
    integer(int32), parameter :: J_r1(N1_terms) = [ &
                                 -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, &
                                 -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, &
                                 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, &
                                 -38, -39, -40, -41]
    real(real64), parameter :: n_r1(N1_terms) = [ &
                               0.14632971213167d0, -0.84548187169114d0, -0.37563603672040d1, &
                               0.33855169168385d1, -0.95791963387872d0, 0.15772038513228d0, &
                               -0.16616417199501d-1, 0.81214629983568d-3, 0.28319080123804d-3, &
                               -0.60706301565874d-3, -0.18990068218419d-1, -0.32529748770505d-1, &
                               -0.21841717175414d-1, -0.52838357969930d-4, -0.47184321073267d-3, &
                               -0.30001780793026d-3, 0.47661393906987d-4, -0.44141845330846d-5, &
                               -0.72694996297594d-15, -0.31679644845054d-4, -0.28270797985312d-5, &
                               -0.85205128120103d-9, -0.22425281908000d-5, -0.65171222895601d-6, &
                               -0.14341729937924d-12, -0.40516996860117d-6, -0.12734301741641d-8, &
                               -0.17424871230634d-9, -0.68762131295531d-18, 0.14478307828521d-19, &
                               0.26335781662795d-22, -0.11947622640071d-22, 0.18228094581404d-23, &
                               -0.93537087292458d-25]

contains
    !> Initialize the IAPWS-97 Region 1 object.
    module pure elemental subroutine initialize_type_iapws97_region1(self)
        implicit none
        !> IAPWS-97 Region 1 objects.
        class(type_iapws97_region1), intent(inout) :: self

        ! Reference values for region 1
        self%T_star = 1386.0d0 ! [K]
        self%p_star = 16.53d6 ! [Pa]
        self%R = 461.526d0 ! [J/(kg·K)]

        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_region1

    !> Calculate the dimensionless Gibbs free energy \(\gamma\) for Region 1.
    module pure elemental subroutine calc_gamma_iapws97_region1(self, tau, pi, property)
        implicit none
        !> IAPWS-97 Region 1 object
        class(type_iapws97_region1), intent(in) :: self
        !> Inverse reduced temperature Tc/T, [-]
        real(real64), intent(in) :: tau
        !> Reduced pressure p/p_c, [-]
        real(real64), intent(in) :: pi
        !> IAPWS Gibbs properties
        type(type_iapws_gamma_property), intent(inout) :: property

        call property%reset()

        ! Calculate gamma and its derivatives
        property%gamma = calc_gamma_region1(tau, pi)
        property%gamma_p = calc_gamma_p_region1(tau, pi)
        property%gamma_t = calc_gamma_t_region1(tau, pi)
        property%gamma_pp = calc_gamma_pp_region1(tau, pi)
        property%gamma_tt = calc_gamma_tt_region1(tau, pi)
        property%gamma_pt = calc_gamma_pt_region1(tau, pi)
    end subroutine calc_gamma_iapws97_region1

    !> Calculate the dimensionless Gibbs free energy \(\gamma\) for Region 1.
    !> Formula: \(\gamma = \sum n_i (7.1 - \pi)^{I_i} (\tau - 1.222)^{J_i}\)
    pure elemental function calc_gamma_region1(tau, pi) result(gamma)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Resulting \(\gamma\) value
        real(real64) :: gamma

        integer(int32) :: i
        real(real64) :: c_gamma

        gamma = 0.0d0
        c_gamma = 0.0d0

        do i = 1, N1_terms
            call kahan_add(gamma, c_gamma, &
                           n_r1(i) * (7.1d0 - pi)**I_r1(i) * (tau - 1.222d0)**J_r1(i))
        end do

    end function calc_gamma_region1

    !> Calculate the first derivative of \(\gamma\) with respect to \(\pi\).
    !> Computes \(\gamma_{\pi} = \left(\frac{\partial \gamma}{\partial \pi}\right)_{\tau}\).
    pure elemental function calc_gamma_p_region1(tau, pi) result(gamma_p)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma_{\pi}\)
        real(real64) :: gamma_p

        integer(int32) :: i
        real(real64) :: term_pi
        real(real64) :: I_val
        real(real64) :: c_gamma_p

        gamma_p = 0.0d0
        c_gamma_p = 0.0d0

        do i = 1, N1_terms
            I_val = real(I_r1(i), real64)
            call kahan_add(gamma_p, c_gamma_p, &
                           -n_r1(i) * I_val * (7.1d0 - pi)**(I_r1(i) - 1) * (tau - 1.222d0)**J_r1(i))
        end do

    end function calc_gamma_p_region1

    !> Calculate the first derivative of \(\gamma\) with respect to \(\tau\).
    !> Computes \(\gamma_{\tau} = \left(\frac{\partial \gamma}{\partial \tau}\right)_{\pi}\).
    pure elemental function calc_gamma_t_region1(tau, pi) result(gamma_t)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma_{\tau}\)
        real(real64) :: gamma_t

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma_t

        gamma_t = 0.0d0
        c_gamma_t = 0.0d0

        do i = 1, N1_terms
            J_val = real(J_r1(i), real64)
            call kahan_add(gamma_t, c_gamma_t, &
                           n_r1(i) * (7.1d0 - pi)**I_r1(i) * J_val * (tau - 1.222d0)**(J_r1(i) - 1))
        end do

    end function calc_gamma_t_region1

    !> Calculate the second derivative of \(\gamma\) with respect to \(\pi\).
    !> Computes \(\gamma_{\pi\pi} = \left(\frac{\partial^2 \gamma}{\partial \pi^2}\right)_{\tau}\).
    pure elemental function calc_gamma_pp_region1(tau, pi) result(gamma_pp)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma_{\pi\pi}\)
        real(real64) :: gamma_pp

        integer(int32) :: i
        real(real64) :: I_val
        real(real64) :: c_gamma_pp

        gamma_pp = 0.0d0
        c_gamma_pp = 0.0d0

        do i = 1, N1_terms
            I_val = real(I_r1(i), real64)
            call kahan_add(gamma_pp, c_gamma_pp, &
                           n_r1(i) * I_val * (I_val - 1.0d0) * (7.1d0 - pi)**(I_r1(i) - 2) * (tau - 1.222d0)**J_r1(i))
        end do

    end function calc_gamma_pp_region1

    !> Calculate the second derivative of \(\gamma\) with respect to \(\tau\).
    !> Computes \(\gamma_{\tau\tau} = \left(\frac{\partial^2 \gamma}{\partial \tau^2}\right)_{\pi}\).
    pure elemental function calc_gamma_tt_region1(tau, pi) result(gamma_tt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma_{\tau\tau}\)
        real(real64) :: gamma_tt

        integer(int32) :: i
        real(real64) :: J_val
        real(real64) :: c_gamma_tt

        gamma_tt = 0.0d0
        c_gamma_tt = 0.0d0

        do i = 1, N1_terms
            J_val = real(J_r1(i), real64)
            call kahan_add(gamma_tt, c_gamma_tt, &
                           n_r1(i) * (7.1d0 - pi)**I_r1(i) * J_val * (J_val - 1.0d0) * (tau - 1.222d0)**(J_r1(i) - 2))
        end do

    end function calc_gamma_tt_region1

    !> Calculate the mixed second derivative of \(\gamma\).
    !> Computes \(\gamma_{\pi\tau} = \frac{\partial^2 \gamma}{\partial \pi \partial \tau}\).
    pure elemental function calc_gamma_pt_region1(tau, pi) result(gamma_pt)
        implicit none
        !> Dimensionless temperature \(\tau\)
        real(real64), intent(in) :: tau
        !> Dimensionless pressure \(\pi\)
        real(real64), intent(in) :: pi
        !> Derivative \(\gamma_{\pi\tau}\)
        real(real64) :: gamma_pt

        integer(int32) :: i
        real(real64) :: I_val, J_val
        real(real64) :: c_gamma_pt

        gamma_pt = 0.0d0
        c_gamma_pt = 0.0d0

        do i = 1, N1_terms
            I_val = real(I_r1(i), real64)
            J_val = real(J_r1(i), real64)
            call kahan_add(gamma_pt, c_gamma_pt, &
                           -n_r1(i) * I_val * (7.1d0 - pi)**(I_r1(i) - 1) * J_val * (tau - 1.222d0)**(J_r1(i) - 1))
        end do
    end function calc_gamma_pt_region1

end submodule iapws97_region1
