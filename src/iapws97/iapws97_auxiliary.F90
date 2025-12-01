submodule(module_iapws97) iapws97_auxiliary
    implicit none
    real(real64), parameter :: nb23(5) = [ &
                               0.34805185628969d3, &
                               -0.11671859879975d1, &
                               0.10192970039326d-2, &
                               0.57254459862746d3, &
                               0.13918839778870d2]
contains
    module pure elemental subroutine initialize_type_iapws97_auxiliary(self)
        implicit none
        class(type_iapws97_auxiliary), intent(inout) :: self

        self%T_star = 1.0d0 ! [K]
        self%p_star = 1.0d6 ! [Pa]

        self%is_initialized = .true.
    end subroutine initialize_type_iapws97_auxiliary

    !> 領域2と3の境界における圧力を計算する（式5）
    !!
    !! IAPWS-IF97 式(5)に基づき、与えられた温度に対する境界圧力を返す。
    !! B23-equation: pi = n1 + n2*theta + n3*theta^2
    !!
    module pure elemental function calc_p_boundary_iapws97_region23(self, temperature) result(pressure)
        implicit none
        !> IAPWS-97 auxiliary model instance
        class(type_iapws97_auxiliary), intent(in) :: self
        !> Temperature [K]
        real(real64), intent(in) :: temperature
        !> Pressure [Pa]
        real(real64) :: pressure

        real(real64) :: theta_val, pi_val

        theta_val = temperature / self%T_star

        pi_val = nb23(1) + theta_val * (nb23(2) + theta_val * nb23(3))

        pressure = pi_val * self%p_star

    end function calc_p_boundary_iapws97_region23

    !> 領域2と3の境界における温度を計算する（式6）
    !!
    !! IAPWS-IF97 式(6)に基づき、与えられた圧力に対する境界温度を返す。
    !! Equation: theta = n4 + sqrt((pi - n5) / n3)
    !!
    module pure elemental function calc_t_boundary_iapws97_region23(self, pressure) result(temperature)
        implicit none
        !> IAPWS-97 auxiliary model instance
        class(type_iapws97_auxiliary), intent(in) :: self
        !> Pressure [Pa]
        real(real64), intent(in) :: pressure
        !> Temperature [K]
        real(real64) :: temperature

        real(real64) :: pi_val, theta_val

        pi_val = pressure / self%p_star

        theta_val = nb23(4) + sqrt(max(0.0d0, (pi_val - nb23(5)) / nb23(3)))

        temperature = theta_val * self%T_star

    end function calc_t_boundary_iapws97_region23
end submodule iapws97_auxiliary
