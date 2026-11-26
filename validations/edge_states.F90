! Emit machine-readable results for the external water comparison.
program edge_states
    use, intrinsic :: iso_fortran_env, only: real64, iostat_end
    use iapws, only: type_water, type_water_state, type_water_saturation, type_water_jacobian, new_water, &
                     WATER_OK, JAC_RHO, JAC_T
    implicit none
    type(type_water) :: water
    type(type_water_state) :: state, inverse
    type(type_water_saturation) :: saturation
    type(type_water_jacobian) :: jacobian
    integer :: formulation, status
    real(real64) :: temperature, pressure
    character(3) :: operation

    do
        read (*,*,iostat=status) formulation, operation, temperature, pressure
        if (status == iostat_end) exit
        if (status /= 0) error stop 'invalid external comparison input'
        water = new_water(formulation)
        select case (operation)
        case ('sat')
            call water%saturation_t(temperature, saturation)
            if (saturation%status /= WATER_OK) then
                write (*,'(i0)') saturation%status
                cycle
            end if
            write (*,'(i0,3(1x,es25.16))') saturation%status, saturation%p, saturation%liquid%rho, saturation%vapor%rho
        case ('tp')
            call water%state_tp_checked(temperature, pressure, state, jac=jacobian)
            if (state%status /= WATER_OK) then
                write (*,'(i0)') state%status
                cycle
            end if
            call water%state_ph_checked(pressure, state%properties%h, inverse)
            write (*,'(i0,6(1x,es25.16),1x,i0,1x,es25.16)') state%status, state%properties%rho, &
                state%properties%cp, state%properties%h, jacobian%ph(JAC_RHO,1), jacobian%ph(JAC_RHO,2), &
                jacobian%ph(JAC_T,2), inverse%status, inverse%properties%T
        case default
            error stop 'invalid external comparison operation'
        end select
    end do
end program edge_states
