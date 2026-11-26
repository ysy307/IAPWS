module module_water_phase_boundaries
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: melting_pressure, sublimation_pressure, ice_ih_melting_slope, sublimation_slope

    integer, parameter, public :: ICE_IH = 4, ICE_III = 5, ICE_V = 6, ICE_VI = 7, ICE_VII = 8
    real(real64), parameter, public :: triple_T = 273.16d0, triple_p = 611.657d0

contains

    ! R14-08(2011), equations (1)-(5). Pressure in Pa.
    pure elemental function melting_pressure(T, solid) result(p)
        implicit none
        real(real64) :: p
        real(real64), intent(in) :: T
        integer, intent(in), optional :: solid
        real(real64) :: theta
        integer :: phase
        p = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        phase = ICE_IH
        if (present(solid)) phase = solid
        select case (phase)
        case (ICE_IH)
            if (T < 251.165d0 .or. T > triple_T) return
            theta = T / triple_T
            p = triple_p * (1d0 + 1195393.37d0 * (1d0 - theta**3) + &
                80818.3159d0 * (1d0 - theta**25.75d0) + 3338.26860d0 * (1d0 - theta**103.75d0))
        case (ICE_III)
            if (T < 251.165d0 .or. T > 256.164d0) return
            p = 208.566d6 * (1d0 + 0.299948d0 * ((T / 251.165d0)**60 - 1d0))
        case (ICE_V)
            if (T < 256.164d0 .or. T > 273.31d0) return
            p = 350.1d6 * (1d0 + 1.18721d0 * ((T / 256.164d0)**8 - 1d0))
        case (ICE_VI)
            if (T < 273.31d0 .or. T > 355d0) return
            p = 632.4d6 * (1d0 + 1.07476d0 * ((T / 273.31d0)**4.6d0 - 1d0))
        case (ICE_VII)
            if (T < 355d0 .or. T > 715d0) return
            theta = T / 355d0
            p = 2216d6 * exp(1.73683d0 * (1d0 - 1d0 / theta) - &
                0.0544606d0 * (1d0 - theta**5) + 8.06106d-8 * (1d0 - theta**22))
        end select
    end function melting_pressure

    ! R14-08(2011), equation (6). Pressure in Pa.
    pure elemental function sublimation_pressure(T) result(p)
        implicit none
        real(real64) :: p
        real(real64), intent(in) :: T
        real(real64) :: theta
        p = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 50d0 .or. T > triple_T) return
        theta = T / triple_T
        p = triple_p * exp((-21.2144006d0 * theta**0.00333333333d0 + &
            27.3203819d0 * theta**1.20666667d0 - 6.10598130d0 * theta**1.70333333d0) / theta)
        if (T == triple_T) p = triple_p
    end function sublimation_pressure

    ! dp/dT on the Ice Ih melting and sublimation curves of R14-08.
    pure elemental function ice_ih_melting_slope(T) result(slope)
        implicit none
        real(real64) :: slope
        real(real64), intent(in) :: T
        real(real64) :: theta
        slope = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 251.165d0 .or. T > triple_T) return
        theta = T / triple_T
        slope = -triple_p / triple_T * (3d0 * 1195393.37d0 * theta**2 + &
            25.75d0 * 80818.3159d0 * theta**24.75d0 + 103.75d0 * 3338.26860d0 * theta**102.75d0)
    end function ice_ih_melting_slope

    pure elemental function sublimation_slope(T) result(slope)
        implicit none
        real(real64) :: slope
        real(real64), intent(in) :: T
        real(real64) :: theta
        slope = ieee_value(0d0, ieee_quiet_nan)
        if (.not. ieee_is_finite(T)) return
        if (T < 50d0 .or. T > triple_T) return
        theta = T / triple_T
        slope = sublimation_pressure(T) / (triple_T * theta**2) * &
            (-21.2144006d0 * (-0.99666666667d0) * theta**0.00333333333d0 + &
              27.3203819d0 * 0.20666667d0 * theta**1.20666667d0 - &
              6.10598130d0 * 0.70333333d0 * theta**1.70333333d0)
    end function sublimation_slope
end module module_water_phase_boundaries
