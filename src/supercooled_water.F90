module module_supercooled_water
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use module_iapws, only: type_iapws_property
    implicit none
    private
    public :: in_supercooled_domain, evaluate_supercooled

    ! IAPWS G12-15, Tables 1 and 2.
    real(real64), parameter :: ci(20) = [ &
        -8.1570681381655d0, 1.2875032d0, 7.0901673598012d0, -3.2779161d-2, 7.3703949d-1, &
        -2.1628622d-1, -5.1782479d0, 4.2293517d-4, 2.3592109d-2, 4.3773754d0, &
        -2.9967770d-3, -9.6558018d-1, 3.7595286d0, 1.2632441d0, 2.8542697d-1, &
        -8.5994947d-1, -3.2916153d-1, 9.0019616d-2, 8.1149726d-2, -3.2788213d0]
    real(real64), parameter :: ai(20) = [ &
        0d0, 0d0, 1d0, -0.2555d0, 1.5762d0, &
        1.6400d0, 3.6385d0, -0.3828d0, 1.6219d0, 4.3287d0, &
        3.4763d0, 5.1556d0, -0.3593d0, 5.0361d0, 2.9786d0, &
        6.2373d0, 4.0460d0, 5.3558d0, 9.0157d0, 1.2194d0]
    real(real64), parameter :: bi(20) = [ &
        0d0, 1d0, 0d0, 2.1051d0, 1.1422d0, &
        0.9510d0, 0d0, 3.6402d0, 2.0760d0, -0.0016d0, &
        2.2769d0, 0.0008d0, 0.3706d0, -0.3975d0, 2.9730d0, &
        -0.3180d0, 2.9805d0, 2.9265d0, 0.4456d0, 0.1298d0]
    real(real64), parameter :: di(20) = [ &
        0d0, 0d0, 0d0, -0.0016d0, 0.6894d0, &
        0.0130d0, 0.0002d0, 0.0435d0, 0.0500d0, 0.0004d0, &
        0.0528d0, 0.0147d0, 0.8584d0, 0.9924d0, 1.0041d0, &
        1.0961d0, 1.0228d0, 1.0303d0, 1.6180d0, 0.5213d0]

contains

    pure elemental function in_supercooled_domain(T, p) result(valid)
        implicit none
        logical :: valid
        real(real64), intent(in) :: T, p
        real(real64) :: q, theta
        valid = .false.
        if (.not. ieee_is_finite(T) .or. .not. ieee_is_finite(p)) return
        if (T <= 0d0 .or. T > 300d0 .or. p < 0d0 .or. p > 400d6) return
        ! G12-15 equations (18) and (19).
        q = p * 1d-6
        if (q < 198.9d0) then
            theta = T / 235.15d0
            valid = q >= 0.1d0 + 228.27d0 * (1d0 - theta**6.243d0) + 15.724d0 * (1d0 - theta**79.81d0)
        else
            valid = T >= 172.82d0 + q * (0.03718d0 + q * (3.403d-5 - 1.573d-8 * q))
        end if
    end function in_supercooled_domain

    pure elemental subroutine evaluate_supercooled(T, p, state, valid)
        implicit none
        real(real64), intent(in) :: T, p
        type(type_iapws_property), intent(inout) :: state
        logical, intent(out) :: valid
        real(real64), parameter :: T0 = 228.2d0, rho0 = 1081.6482d0, R = 461.523087d0, &
            G0 = R * T0, p0 = rho0 * G0, omega0 = 0.5212269d0, L0 = 0.76317954d0, &
            k0 = 0.072158686d0, k1 = -0.31569232d0, k2 = 5.2992608d0, K = sqrt(1d0 + k2**2)
        real(real64) :: y, z, q, tau, a, b, c, d, term, bt, bp, btt, btp, bpp, background
        real(real64) :: root, omega, L, Lt, Lp, Ltt, Ltp, Lpp, x, lo, hi, f, df, trial, Fmix, Fp
        real(real64) :: g, gt, gp, gtt, gpp, gtp, asym
        integer :: i

        valid = .false.
        y = T / T0
        tau = y - 1d0
        q = p / p0
        z = q + 300d6 / p0
        asym = k1 * (q - k2 * tau)
        root = sqrt((1d0 + k0 * k2 + asym)**2 - 4d0 * k0 * k2 * asym)
        L = L0 * K * (1d0 + k0 * k2 + k1 * (q + k2 * tau) - root) / (2d0 * k1 * k2)
        omega = 2d0 + omega0 * q

        ! Brackets from Table 4 and safeguarded Newton on equation (8).
        if (omega < (10d0 / 9d0) * (log(19d0) - L)) then
            lo = 0.049d0
            hi = 0.5d0
        else if (omega < (50d0 / 49d0) * (log(99d0) - L)) then
            lo = 0.0099d0
            hi = 0.051d0
        else
            lo = 0.99d0 * exp(-(50d0 / 49d0) * L - omega)
            hi = min(0.0101d0, 1.1d0 * exp(-L - omega))
        end if
        x = (lo + hi) * 0.5d0
        do i = 1, 60
            f = L + log(x / (1d0 - x)) + omega * (1d0 - 2d0 * x)
            df = 1d0 / (x * (1d0 - x)) - 2d0 * omega
            if (abs(f) <= 8d0 * epsilon(x)) exit
            if (f < 0d0) then
                lo = x
            else
                hi = x
            end if
            trial = (lo + hi) * 0.5d0
            if (df > 0d0) then
                if (x - f / df > lo .and. x - f / df < hi) trial = x - f / df
            end if
            if (trial == x) exit
            x = trial
        end do
        if (abs(f) > 32d0 * epsilon(x) .or. df <= 0d0) return

        background = 0d0
        bt = 0d0
        bp = 0d0
        btt = 0d0
        btp = 0d0
        bpp = 0d0
        do i = 1, size(ci)
            a = ai(i)
            b = bi(i)
            c = ci(i)
            d = di(i)
            term = c * y**a * z**b * exp(-d * z)
            background = background + term
            bt = bt + term * a
            bp = bp + term * (b - d * z)
            btt = btt + term * a * (a - 1d0)
            btp = btp + term * a * (b - d * z)
            bpp = bpp + term * ((b - d * z)**2 - b)
        end do
        Lt = L0 * K * (1d0 + (1d0 - k0 * k2 + asym) / root) * 0.5d0
        Lp = L0 * K * (root + k0 * k2 - asym - 1d0) / (2d0 * k2 * root)
        Lpp = -2d0 * L0 * K * k0 * k1 / root**3
        Ltp = -k2 * Lpp
        Ltt = k2**2 * Lpp
        Fmix = x * L + x * log(x) + (1d0 - x) * log(1d0 - x) + omega * x * (1d0 - x)
        Fp = x * Lp + omega0 * x * (1d0 - x)
        asym = Lp + omega0 * (1d0 - 2d0 * x)
        g = G0 * (background + y * Fmix)
        gt = R * (bt / y + Fmix + y * x * Lt)
        gp = (bp / z + y * Fp) / rho0
        gtt = R / T0 * (btt / y**2 + 2d0 * x * Lt + y * (x * Ltt - Lt**2 / df))
        gtp = (btp / (y * z) + Fp + y * (x * Ltp - Lt * asym / df)) / (rho0 * T0)
        gpp = (bpp / z**2 + y * (x * Lpp - asym**2 / df)) / (rho0 * p0)
        if (gp <= 0d0 .or. gtt >= 0d0 .or. gpp >= 0d0) return

        state%T = T
        state%p = p
        state%nu = gp
        state%rho = 1d0 / gp
        state%s = -gt
        state%h = g - T * gt
        state%u = state%h - p * gp
        state%cp = -T * gtt
        state%cv = state%cp + T * gtp**2 / gpp
        state%alpha = gtp / gp
        state%beta = -gtp / gpp
        state%kappa_T = -gpp / gp
        state%kappa_s = state%kappa_T * (state%cv / state%cp)
        state%drho_dT = -gtp / gp**2
        state%drho_dP = -gpp / gp**2
        if (state%kappa_s <= 0d0) return
        state%w = sqrt(gp / state%kappa_s)
        valid = .true.
    end subroutine evaluate_supercooled
end module module_supercooled_water
