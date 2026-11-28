submodule(module_iapws97) iapws97_saturation_batch
    use module_iapws_batch, only: IAPWS_BATCH_OK, IAPWS_BATCH_INVALID, IAPWS_BATCH_TILE
    use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
contains
    module pure recursive subroutine calc_saturation_properties_batch_iapws97(self, T_in, rho_vap, rho_liq, &
                                                                              cp_vap, cp_liq, drho_dT_vap, &
                                                                              drho_dT_liq, drho_dP_vap, &
                                                                              drho_dP_liq, density_dT_vap, &
                                                                              density_dT_liq, dcp_dT_vap, &
                                                                              dcp_dT_liq, latent_heat, dlatent_heat_dT, status)
        implicit none
        class(type_iapws97), intent(in) :: self
        real(real64), intent(in) :: T_in(:)
        real(real64), intent(out), optional :: rho_vap(:), rho_liq(:), cp_vap(:), cp_liq(:), drho_dT_vap(:), &
                                               drho_dT_liq(:), drho_dP_vap(:), drho_dP_liq(:), density_dT_vap(:), &
                                               density_dT_liq(:), &
                                               dcp_dT_vap(:), dcp_dT_liq(:), latent_heat(:), dlatent_heat_dT(:)
        integer, intent(out), optional :: status(:)
        integer :: indices(IAPWS_BATCH_TILE), count, j
        real(real64) :: nan, values(1, 14)
        logical :: valid, want_values, scalar_model
        if (present(rho_vap)) then
            if (size(rho_vap) /= size(T_in)) error stop 'IAPWS batch: saturation rho_vap size'
        end if
        if (present(rho_liq)) then
            if (size(rho_liq) /= size(T_in)) error stop 'IAPWS batch: saturation rho_liq size'
        end if
        if (present(cp_vap)) then
            if (size(cp_vap) /= size(T_in)) error stop 'IAPWS batch: saturation cp_vap size'
        end if
        if (present(cp_liq)) then
            if (size(cp_liq) /= size(T_in)) error stop 'IAPWS batch: saturation cp_liq size'
        end if
        if (present(drho_dT_vap)) then
            if (size(drho_dT_vap) /= size(T_in)) error stop 'IAPWS batch: saturation drho_dT_vap size'
        end if
        if (present(drho_dT_liq)) then
            if (size(drho_dT_liq) /= size(T_in)) error stop 'IAPWS batch: saturation drho_dT_liq size'
        end if
        if (present(drho_dP_vap)) then
            if (size(drho_dP_vap) /= size(T_in)) error stop 'IAPWS batch: saturation drho_dP_vap size'
        end if
        if (present(drho_dP_liq)) then
            if (size(drho_dP_liq) /= size(T_in)) error stop 'IAPWS batch: saturation drho_dP_liq size'
        end if
        if (present(density_dT_vap)) then
            if (size(density_dT_vap) /= size(T_in)) error stop 'IAPWS batch: saturation density_dT_vap size'
        end if
        if (present(density_dT_liq)) then
            if (size(density_dT_liq) /= size(T_in)) error stop 'IAPWS batch: saturation density_dT_liq size'
        end if
        if (present(dcp_dT_vap)) then
            if (size(dcp_dT_vap) /= size(T_in)) error stop 'IAPWS batch: saturation dcp_dT_vap size'
        end if
        if (present(dcp_dT_liq)) then
            if (size(dcp_dT_liq) /= size(T_in)) error stop 'IAPWS batch: saturation dcp_dT_liq size'
        end if
        if (present(latent_heat)) then
            if (size(latent_heat) /= size(T_in)) error stop 'IAPWS batch: saturation latent_heat size'
        end if
        if (present(dlatent_heat_dT)) then
            if (size(dlatent_heat_dT) /= size(T_in)) error stop 'IAPWS batch: saturation dlatent_heat_dT size'
        end if
        if (present(status)) then
            if (size(status) /= size(T_in)) error stop 'IAPWS batch: saturation status size'
        end if
        want_values = any([present(rho_vap), present(rho_liq), present(cp_vap), present(cp_liq), &
                           present(drho_dT_vap), present(drho_dT_liq), present(drho_dP_vap), present(drho_dP_liq), &
                           present(density_dT_vap), present(density_dT_liq), present(dcp_dT_vap), present(dcp_dT_liq), &
                           present(latent_heat), present(dlatent_heat_dT)])
        select type (self)
        type is (type_iapws97)
            scalar_model = .false.
        class default
            scalar_model = .true.
        end select
        scalar_model = scalar_model .or. size(T_in) < 4
        nan = ieee_value(0d0, ieee_quiet_nan); count = 0
        do j = 1, size(T_in)
            valid = ieee_is_finite(T_in(j))
            if (valid) valid = T_in(j) >= 273.15d0 .and. T_in(j) <= 647.096d0
            if (.not. valid) then
                if (present(status)) status(j) = IAPWS_BATCH_INVALID
                if (present(rho_vap)) rho_vap(j) = nan
                if (present(rho_liq)) rho_liq(j) = nan
                if (present(cp_vap)) cp_vap(j) = nan
                if (present(cp_liq)) cp_liq(j) = nan
                if (present(drho_dT_vap)) drho_dT_vap(j) = nan
                if (present(drho_dT_liq)) drho_dT_liq(j) = nan
                if (present(drho_dP_vap)) drho_dP_vap(j) = nan
                if (present(drho_dP_liq)) drho_dP_liq(j) = nan
                if (present(density_dT_vap)) density_dT_vap(j) = nan
                if (present(density_dT_liq)) density_dT_liq(j) = nan
                if (present(dcp_dT_vap)) dcp_dT_vap(j) = nan
                if (present(dcp_dT_liq)) dcp_dT_liq(j) = nan
                if (present(latent_heat)) latent_heat(j) = nan
                if (present(dlatent_heat_dT)) dlatent_heat_dT(j) = nan
                cycle
            end if
            if (present(status)) status(j) = IAPWS_BATCH_OK
            if (.not. want_values .and. .not. (present(status) .and. T_in(j) > 623.15d0)) cycle
            if (scalar_model .or. T_in(j) > 623.15d0) then
                ! Keep critical values, failed-pair NaNs, and independent
                ! region-3 root histories exactly as the scalar evaluator.
                if (present(dcp_dT_vap) .and. present(dcp_dT_liq)) then
                    call self%calc_saturation_properties(T_in(j), rho_vap=values(1, 1), rho_liq=values(1, 2), &
                                                         cp_vap=values(1, 3), cp_liq=values(1, 4), drho_dT_vap=values(1, 5), &
                                                         drho_dT_liq=values(1, 6), drho_dP_vap=values(1, 7), &
                                                         drho_dP_liq=values(1, 8), density_dT_vap=values(1, 9), &
                                                         density_dT_liq=values(1, 10), dcp_dT_vap=values(1, 11), &
                                                         dcp_dT_liq=values(1, 12), latent_heat=values(1, 13), &
                                                         dlatent_heat_dT=values(1, 14))
                else if (present(dcp_dT_vap)) then
                    call self%calc_saturation_properties(T_in(j), rho_vap=values(1, 1), rho_liq=values(1, 2), &
                                                         cp_vap=values(1, 3), cp_liq=values(1, 4), drho_dT_vap=values(1, 5), drho_dT_liq=values(1, 6), &
                                                         drho_dP_vap=values(1, 7), drho_dP_liq=values(1, 8), density_dT_vap=values(1, 9), &
                                                         density_dT_liq=values(1, 10), dcp_dT_vap=values(1, 11), latent_heat=values(1, 13), &
                                                         dlatent_heat_dT=values(1, 14))
                else if (present(dcp_dT_liq)) then
                    call self%calc_saturation_properties(T_in(j), rho_vap=values(1, 1), rho_liq=values(1, 2), &
                                                         cp_vap=values(1, 3), cp_liq=values(1, 4), drho_dT_vap=values(1, 5), drho_dT_liq=values(1, 6), &
                                                         drho_dP_vap=values(1, 7), drho_dP_liq=values(1, 8), density_dT_vap=values(1, 9), &
                                                         density_dT_liq=values(1, 10), dcp_dT_liq=values(1, 12), latent_heat=values(1, 13), &
                                                         dlatent_heat_dT=values(1, 14))
                else
                    call self%calc_saturation_properties(T_in(j), rho_vap=values(1, 1), rho_liq=values(1, 2), &
                                                         cp_vap=values(1, 3), cp_liq=values(1, 4), drho_dT_vap=values(1, 5), drho_dT_liq=values(1, 6), &
                                                         drho_dP_vap=values(1, 7), drho_dP_liq=values(1, 8), density_dT_vap=values(1, 9), &
                                                         density_dT_liq=values(1, 10), latent_heat=values(1, 13), dlatent_heat_dT=values(1, 14))
                end if
                if (present(status)) then
                    if (.not. all(ieee_is_finite(values(1, :2)))) status(j) = IAPWS_BATCH_INVALID
                end if
                call store_saturation([j], values, rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, drho_dT_liq, &
                                      drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, latent_heat, &
                                      dlatent_heat_dT)
                cycle
            end if
            count = count + 1; indices(count) = j
            if (count == IAPWS_BATCH_TILE) then
                call saturation_low(self, indices, T_in, rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, &
                                    drho_dT_liq, drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, &
                                    latent_heat, dlatent_heat_dT)
                count = 0
            end if
        end do
        if (count > 0) call saturation_low(self, indices(:count), T_in, rho_vap, rho_liq, cp_vap, cp_liq, &
                                           drho_dT_vap, drho_dT_liq, drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, &
                                           dcp_dT_liq, latent_heat, dlatent_heat_dT)
    end subroutine

    pure recursive subroutine saturation_low(self, indices, T_in, rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, &
                                             drho_dT_liq, drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, &
                                             latent_heat, dlatent_heat_dT)
        implicit none
        class(type_iapws97), intent(in) :: self
        integer, intent(in) :: indices(:)
        real(real64), intent(in) :: T_in(:)
        real(real64), intent(inout), optional :: rho_vap(:), rho_liq(:), cp_vap(:), cp_liq(:), drho_dT_vap(:), &
                                                 drho_dT_liq(:), drho_dP_vap(:), drho_dP_liq(:), density_dT_vap(:), density_dT_liq(:), &
                                                 dcp_dT_vap(:), dcp_dT_liq(:), latent_heat(:), dlatent_heat_dT(:)
        real(real64) :: t(IAPWS_BATCH_TILE), pressure(IAPWS_BATCH_TILE), slope(IAPWS_BATCH_TILE)
        real(real64) :: gv(IAPWS_BATCH_TILE, 6), gl(IAPWS_BATCH_TILE, 6), values(IAPWS_BATCH_TILE, 14)
        real(real64) :: hv, hl, hpv, hpl, ct, cp_pressure
        integer :: j, n
        n = size(indices)
        do j = 1, n
            t(j) = T_in(indices(j))
            pressure(j) = self%region4%calc_psat(t(j)); slope(j) = self%region4%calc_dpsat_dT(t(j))
        end do
        call self%region2%calc_gamma_batch(t(:n), pressure(:n), gv(:n, :))
        call self%region1%calc_gamma_batch(t(:n), pressure(:n), gl(:n, :))
#if !defined(__NVCOMPILER)
        !$omp simd private(hv,hl,hpv,hpl)
#endif
        do j = 1, n
            values(j, 1) = 1d0 / gv(j, 2); values(j, 2) = 1d0 / gl(j, 2)
            values(j, 3) = -t(j) * gv(j, 5); values(j, 4) = -t(j) * gl(j, 5)
            values(j, 5) = -gv(j, 6) / gv(j, 2)**2; values(j, 6) = -gl(j, 6) / gl(j, 2)**2
            values(j, 7) = -gv(j, 4) / gv(j, 2)**2; values(j, 8) = -gl(j, 4) / gl(j, 2)**2
            values(j, 9) = values(j, 5) + values(j, 7) * slope(j)
            values(j, 10) = values(j, 6) + values(j, 8) * slope(j)
            hv = gv(j, 1) - t(j) * gv(j, 3); hl = gl(j, 1) - t(j) * gl(j, 3)
            hpv = gv(j, 2) - t(j) * gv(j, 6); hpl = gl(j, 2) - t(j) * gl(j, 6)
            values(j, 13) = hv - hl
            values(j, 14) = (values(j, 3) - values(j, 4)) + (hpv - hpl) * slope(j)
        end do
        if (present(dcp_dT_vap)) then
            do j = 1, n
                call self%region2%calc_cp_derivatives(t(j), pressure(j), ct, cp_pressure)
                values(j, 11) = ct + cp_pressure * slope(j)
            end do
        end if
        if (present(dcp_dT_liq)) then
            do j = 1, n
                call self%region1%calc_cp_derivatives(t(j), pressure(j), ct, cp_pressure)
                values(j, 12) = ct + cp_pressure * slope(j)
            end do
        end if
        call store_saturation(indices, values(:n, :), rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, drho_dT_liq, &
                              drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, latent_heat, &
                              dlatent_heat_dT)
    end subroutine

    pure recursive subroutine store_saturation(indices, values, rho_vap, rho_liq, cp_vap, cp_liq, drho_dT_vap, &
                                               drho_dT_liq, drho_dP_vap, drho_dP_liq, density_dT_vap, density_dT_liq, dcp_dT_vap, dcp_dT_liq, &
                                               latent_heat, dlatent_heat_dT)
        implicit none
        integer, intent(in) :: indices(:)
        real(real64), intent(in) :: values(:, :)
        real(real64), intent(inout), optional :: rho_vap(:), rho_liq(:), cp_vap(:), cp_liq(:), drho_dT_vap(:), &
                                                 drho_dT_liq(:), drho_dP_vap(:), drho_dP_liq(:), density_dT_vap(:), density_dT_liq(:), &
                                                 dcp_dT_vap(:), dcp_dT_liq(:), latent_heat(:), dlatent_heat_dT(:)
        integer :: j
        if (present(rho_vap)) then
            do j = 1, size(indices)
                rho_vap(indices(j)) = values(j, 1)
            end do
        end if
        if (present(rho_liq)) then
            do j = 1, size(indices)
                rho_liq(indices(j)) = values(j, 2)
            end do
        end if
        if (present(cp_vap)) then
            do j = 1, size(indices)
                cp_vap(indices(j)) = values(j, 3)
            end do
        end if
        if (present(cp_liq)) then
            do j = 1, size(indices)
                cp_liq(indices(j)) = values(j, 4)
            end do
        end if
        if (present(drho_dT_vap)) then
            do j = 1, size(indices)
                drho_dT_vap(indices(j)) = values(j, 5)
            end do
        end if
        if (present(drho_dT_liq)) then
            do j = 1, size(indices)
                drho_dT_liq(indices(j)) = values(j, 6)
            end do
        end if
        if (present(drho_dP_vap)) then
            do j = 1, size(indices)
                drho_dP_vap(indices(j)) = values(j, 7)
            end do
        end if
        if (present(drho_dP_liq)) then
            do j = 1, size(indices)
                drho_dP_liq(indices(j)) = values(j, 8)
            end do
        end if
        if (present(density_dT_vap)) then
            do j = 1, size(indices)
                density_dT_vap(indices(j)) = values(j, 9)
            end do
        end if
        if (present(density_dT_liq)) then
            do j = 1, size(indices)
                density_dT_liq(indices(j)) = values(j, 10)
            end do
        end if
        if (present(dcp_dT_vap)) then
            do j = 1, size(indices)
                dcp_dT_vap(indices(j)) = values(j, 11)
            end do
        end if
        if (present(dcp_dT_liq)) then
            do j = 1, size(indices)
                dcp_dT_liq(indices(j)) = values(j, 12)
            end do
        end if
        if (present(latent_heat)) then
            do j = 1, size(indices)
                latent_heat(indices(j)) = values(j, 13)
            end do
        end if
        if (present(dlatent_heat_dT)) then
            do j = 1, size(indices)
                dlatent_heat_dT(indices(j)) = values(j, 14)
            end do
        end if
    end subroutine
end submodule iapws97_saturation_batch
