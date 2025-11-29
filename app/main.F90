program main
    use, intrinsic :: iso_fortran_env
    use :: iapws95_interface
    implicit none

    real(real64) :: input_value, output_value
    input_value = 3.0d0
    output_value = example_function(input_value)

    print *, "Input Value: ", input_value
    print *, "Output Value: ", output_value

end program main
