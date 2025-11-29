module iapws95_interface
    use, intrinsic :: iso_fortran_env
    implicit none

contains
    function example_function(x) result(y)
        real(real64), intent(in) :: x
        real(real64) :: y

        ! Example implementation: y = x^2
        y = x * x
    end function example_function

end module iapws95_interface
