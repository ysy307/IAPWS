module module_water_transport
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: viscosity, thermal_conductivity, surface_tension, industrial_reference_slope

    real(real64), parameter :: Tc = 647.096d0, rhoc = 322d0, pc = 22.064d6, &
        pi = 3.14159265358979323846d0
    ! R12-08 Table 2 and R15-11 Table 2, grouped for nested Horner evaluation.
    real(real64), parameter :: H(0:6,0:5) = reshape([ &
        0.520094d0, 0.222531d0, -0.281378d0, 0.161913d0, -0.0325372d0, 0.0d0, 0.0d0, &
        0.0850895d0, 0.999115d0, -0.906851d0, 0.257399d0, 0.0d0, 0.0d0, 0.0d0, &
        -1.08374d0, 1.88797d0, -0.772479d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        -0.289555d0, 1.26613d0, -0.489837d0, 0.0d0, 0.0698452d0, 0.0d0, -0.00435673d0, &
        0.0d0, 0.0d0, -0.25704d0, 0.0d0, 0.0d0, 0.00872102d0, 0.0d0, &
        0.0d0, 0.120573d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, -0.000593264d0], [7,6])
    real(real64), parameter :: L(0:5,0:4) = reshape([ &
        1.60397357d0, -0.646013523d0, 0.111443906d0, 0.102997357d0, -0.0504123634d0, 0.00609859258d0, &
        2.33771842d0, -2.78843778d0, 1.53616167d0, -0.463045512d0, 0.0832827019d0, -0.00719201245d0, &
        2.19650529d0, -4.54580785d0, 3.55777244d0, -1.40944978d0, 0.275418278d0, -0.0205938816d0, &
        -1.21051378d0, 1.60812989d0, -0.621178141d0, 0.0716373224d0, 0.0d0, 0.0d0, &
        -2.720337d0, 4.57586331d0, -3.18369245d0, 1.1168348d0, -0.19268305d0, 0.012913842d0], [6,5])

contains

    pure elemental function viscosity(T, rho, drho_dp, reference_drho_dp) result(mu)
        implicit none
        real(real64) :: mu
        real(real64), intent(in) :: T, rho
        real(real64), intent(in), optional :: drho_dp, reference_drho_dp
        real(real64) :: tau, d, delta, x, qcx, qdx, angle, w, lw, y

        mu = ieee_value(0d0, ieee_quiet_nan)
        if (.not. valid_input(T, rho)) return
        tau = Tc / T
        d = rho / rhoc
        delta = d - 1d0
        ! R12-08 equations (11), (12), and (14)-(21).
        mu = 1d-4 * exp(d * evaluate_polynomial(H, delta, tau - 1d0)) / &
            (sqrt(tau) * (1.67752d0 + tau * (2.20462d0 + tau * (0.6366564d0 - 0.241605d0 * tau))))
        if (.not. present(drho_dp) .or. .not. present(reference_drho_dp)) return
        x = correlation_length(d, tau, drho_dp, reference_drho_dp)
        qcx = x / 1.9d0
        qdx = x / 1.1d0
        if (x <= 0.3817016416d0) then
            y = qcx / 5d0 * qdx**5 * (1d0 - qcx + qcx**2 - (765d0 / 504d0) * qdx**2)
        else
            angle = atan(qdx)
            w = sqrt(abs((qcx - 1d0) / (qcx + 1d0))) * tan(angle * 0.5d0)
            lw = 2d0 * atan(w)
            if (qcx > 1d0) lw = log((1d0 + w) / (1d0 - w))
            y = sin(3d0 * angle) / 12d0 - sin(2d0 * angle) / (4d0 * qcx) + &
                (1d0 / qcx**2 - 1.25d0) * sin(angle) - &
                ((1d0 - 1.5d0 * qcx**2) * angle - abs(qcx**2 - 1d0)**1.5d0 * lw) / qcx**3
        end if
        mu = mu * exp(0.068d0 * y)
    end function viscosity

    pure elemental function thermal_conductivity(T, rho, cp, cv, mu, drho_dp, reference_drho_dp) result(k)
        implicit none
        real(real64) :: k
        real(real64), intent(in) :: T, rho
        real(real64), intent(in), optional :: cp, cv, mu, drho_dp, reference_drho_dp
        real(real64) :: tau, d, y, q, z, ratio, ref

        k = ieee_value(0d0, ieee_quiet_nan)
        if (.not. valid_input(T, rho)) return
        tau = Tc / T
        d = rho / rhoc
        ! R15-11 equations (16)-(22), SI units.
        k = 1d-3 * exp(d * evaluate_polynomial(L, d - 1d0, tau - 1d0)) / &
            (sqrt(tau) * (2.443221d-3 + tau * (1.323095d-2 + tau * &
            (6.770357d-3 + tau * (-3.454586d-3 + 4.096266d-4 * tau)))))
        if (.not. present(cp) .or. .not. present(cv) .or. .not. present(mu) .or. .not. present(drho_dp)) return
        if (d <= 0d0) return
        if (.not. ieee_is_finite(cp) .or. .not. ieee_is_finite(cv) .or. .not. ieee_is_finite(mu)) return
        if (cp <= 0d0 .or. cv <= 0d0 .or. mu <= 0d0) return
        if (present(reference_drho_dp)) then
            ref = reference_drho_dp
        else
            ref = industrial_reference_slope(rho)
        end if
        y = correlation_length(d, tau, drho_dp, ref) / 0.4d0
        if (y < 1.2d-7) return
        ratio = cv / cp
        q = y / (1d0 + y**3 / (3d0 * d**2))
        ! Evaluate 1-exp(-q) without cancellation for small q.
        z = 1d0 - exp(-q)
        if (q < 1d-4) z = q * (1d0 + q * (-0.5d0 + q * (1d0 / 6d0 - q / 24d0)))
        z = (2d0 / (pi * y)) * ((1d0 - ratio) * atan(y) + y * ratio - z)
        k = k + 177.8514d-9 * d * cp / 461.51805d0 / (tau * mu) * z
    end function thermal_conductivity

    pure elemental function industrial_reference_slope(rho) result(slope)
        implicit none
        real(real64) :: slope
        real(real64), intent(in) :: rho
        real(real64), parameter :: a(0:5,5) = reshape([ &
            6.53786807199516d0, -5.61149954923348d0, 3.39624167361325d0, &
            -2.27492629730878d0, 10.2631854662709d0, 1.97815050331519d0, &
            6.52717759281799d0, -6.30816983387575d0, 8.08379285492595d0, &
            -9.82240510197603d0, 12.1358413791395d0, -5.54349664571295d0, &
            5.35500529896124d0, -3.96415689925446d0, 8.91990208918795d0, &
            -12.0338729505790d0, 9.19494865194302d0, -2.16866274479712d0, &
            1.55225959906681d0, 0.464621290821181d0, 8.93237374861479d0, &
            -11.0321960061126d0, 6.16780999933360d0, -0.965458722086812d0, &
            1.11999926419994d0, 0.595748562571649d0, 9.88952565078920d0, &
            -10.3255051147040d0, 4.66861294457414d0, -0.503243546373828d0], [6,5])
        real(real64) :: d, polynomial_sum
        integer :: i, j
        d = rho / rhoc
        j = 1 + count(d > [0.310559006d0, 0.776397516d0, 1.242236025d0, 1.863354037d0])
        polynomial_sum = a(5,j)
        do i = 4, 0, -1
            polynomial_sum = a(i,j) + d * polynomial_sum
        end do
        slope = (rhoc / pc) / polynomial_sum
    end function industrial_reference_slope

    pure elemental function correlation_length(d, tau, slope, reference) result(x)
        implicit none
        real(real64) :: x
        real(real64), intent(in) :: d, tau, slope, reference
        real(real64) :: enhancement
        x = 0d0
        if (.not. ieee_is_finite(slope) .or. .not. ieee_is_finite(reference)) return
        enhancement = (pc / rhoc) * d * (slope - 1.5d0 * tau * reference)
        if (enhancement > 0d0) x = 0.13d0 * (enhancement / 0.06d0)**(0.63d0 / 1.239d0)
    end function correlation_length

    pure function evaluate_polynomial(coefficients, x, y) result(polynomial_sum)
        implicit none
        real(real64) :: polynomial_sum
        real(real64), intent(in) :: coefficients(0:,0:), x, y
        real(real64) :: row
        integer :: i, j
        polynomial_sum = 0d0
        do j = ubound(coefficients, 2), 0, -1
            row = coefficients(ubound(coefficients, 1), j)
            do i = ubound(coefficients, 1) - 1, 0, -1
                row = coefficients(i,j) + x * row
            end do
            polynomial_sum = row + y * polynomial_sum
        end do
    end function evaluate_polynomial

    pure elemental function valid_input(T, rho) result(is_valid)
        implicit none
        logical :: is_valid
        real(real64), intent(in) :: T, rho
        is_valid = .false.
        if (.not. ieee_is_finite(T) .or. .not. ieee_is_finite(rho)) return
        is_valid = T >= 248.15d0 .and. T <= 1173.15d0 .and. rho >= 0d0 .and. rho <= 1800d0
    end function valid_input

    pure elemental function surface_tension(T) result(sigma)
        implicit none
        real(real64) :: sigma
        real(real64), intent(in) :: T
        real(real64) :: theta
        sigma = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 248.15d0 .or. T > Tc) return
        theta = 1d0 - T / Tc
        sigma = 0.2358d0 * theta**1.256d0 * (1d0 - 0.625d0 * theta)
    end function surface_tension
end module module_water_transport
