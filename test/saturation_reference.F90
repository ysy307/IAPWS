! Independent fixtures from jjgomera/iapws 1.5.5 (https://github.com/jjgomera/iapws).
! Generated with _PSat_T, _Region1, _Region2 and _Region3, SI conversions applied.
! Above 623.15 K, SciPy brentq refines both stable EOS density roots (xtol=1e-10);
! the package's unrefined backward volume approximation is not used.
! Partial density derivatives use -rho*alfav and rho*kt/1e6 respectively.
! These fixed fixtures require no Python, SciPy, or other runtime dependency.
program saturation_reference
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use iapws, only: type_iapws97
    implicit none
    integer, parameter :: n=8
    real(real64), parameter :: temperatures(n)=[300d0,500d0,623.15d0,623.1501d0,630d0,640d0,646d0,647.09d0]
    real(real64), parameter :: rho_vap(n) = [ &
        2.5587188671471278e-02_real64, 1.3197636894926514e+01_real64, &
        1.1362433067027342e+02_real64, 1.1361295025469194e+02_real64, &
        1.3289447773989443e+02_real64, 1.7740124274998382e+02_real64, &
        2.4464342184993683e+02_real64, 3.1484898416451972e+02_real64 ]
    real(real64), parameter :: cp_vap(n) = [ &
        1.9139326762186347e+03_real64, 3.4626212454950955e+03_real64, &
        1.6641479571884935e+04_real64, 1.6692305150736964e+04_real64, &
        2.2788506380717405e+04_real64, 5.2520818475767366e+04_real64, &
        3.4402414332112507e+05_real64, 7.1255648455855921e+07_real64 ]
    real(real64), parameter :: rho_T_vap(n) = [ &
        -8.6393077095533162e-05_real64, -4.9852028215928670e-02_real64, &
        -2.5746787894456542e+00_real64, -2.5801555805869407e+00_real64, &
        -4.2779055614845207e+00_real64, -1.4208398504327413e+01_real64, &
        -1.4170677131656990e+02_real64, -4.1308777252526088e+04_real64 ]
    real(real64), parameter :: rho_p_vap(n) = [ &
        7.2484586293153850e-06_real64, 5.9883579675213257e-06_real64, &
        2.4772737774646632e-05_real64, 2.4805621764072873e-05_real64, &
        3.4576355603710482e-05_real64, 8.5087906312676783e-05_real64, &
        6.3859680469053009e-04_real64, 1.5633396386524326e-01_real64 ]
    real(real64), parameter :: rho_liq(n) = [ &
        9.9651426293008910e+02_real64, 8.3131795923703100e+02_real64, &
        5.7468934171386263e+02_real64, 5.7466999666558581e+02_real64, &
        5.4432837706147416e+02_real64, 4.8161217221244669e+02_real64, &
        4.0025730181175373e+02_real64, 3.2927680215661229e+02_real64 ]
    real(real64), parameter :: cp_liq(n) = [ &
        4.1813730858029867e+03_real64, 4.6590182583982851e+03_real64, &
        1.0101976664061616e+04_real64, 1.0121686720129210e+04_real64, &
        1.2762541834471154e+04_real64, 2.6468397119692880e+04_real64, &
        1.9529583324543817e+05_real64, 6.0822996465536430e+07_real64 ]
    real(real64), parameter :: rho_T_liq(n) = [ &
        -2.7331805638130080e-01_real64, -1.3679570737758737e+00_real64, &
        -5.9565660919016636e+00_real64, -5.9768588682115853e+00_real64, &
        -8.1727820347800488e+00_real64, -1.9432856975821345e+01_real64, &
        -1.4618895208177494e+02_real64, -3.7408309457584226e+04_real64 ]
    real(real64), parameter :: rho_p_liq(n) = [ &
        4.4808403089087112e-07_real64, 9.4193485345001362e-07_real64, &
        9.6314516542162968e-06_real64, 9.6727705608552609e-06_real64, &
        1.4906439996302911e-05_real64, 4.5334644470503741e-05_real64, &
        4.5049535473197576e-04_real64, 1.3732288271123833e-01_real64 ]
    real(real64), parameter :: latent(n) = [ &
        2.4373180174949192e+06_real64, 1.8271251138824497e+06_real64, &
        8.9273378561382857e+05_real64, 8.9274032549461804e+05_real64, &
        7.8009052770228707e+05_real64, 5.5243239819652215e+05_real64, &
        2.6868329271274753e+05_real64, 2.4155375096132047e+04_real64 ]
    type(type_iapws97) :: water
    real(real64) :: vapor(n), liquid(n), values(n)
    call water%initialize()
    call water%calc_saturation_density(temperatures,vapor,liquid)
    call check(vapor,rho_vap,'vapor density')
    call check(liquid,rho_liq,'liquid density')
    if (any(vapor>=liquid)) error stop 'saturation phases must be distinct and ordered'
    call water%calc_saturation_cp(temperatures,vapor,liquid)
    call check(vapor,cp_vap,'vapor cp')
    call check(liquid,cp_liq,'liquid cp')
    call water%calc_latent_heat(values,T_in=temperatures)
    call check(values,latent,'latent heat')
    ! These are the existing partial derivatives evaluated at saturation states,
    ! rather than derivatives following the saturation pressure curve.
    call water%calc_saturation_drho_dT(temperatures,vapor,liquid)
    call check(vapor,rho_T_vap,'vapor density temperature partial')
    call check(liquid,rho_T_liq,'liquid density temperature partial')
    call water%calc_saturation_drho_dP(temperatures,vapor,liquid)
    call check(vapor,rho_p_vap,'vapor density pressure partial')
    call check(liquid,rho_p_liq,'liquid density pressure partial')
    call water%calc_saturation_properties(temperatures,rho_vap=vapor,rho_liq=liquid)
    call check(vapor,rho_vap,'bulk vapor density')
    call check(liquid,rho_liq,'bulk liquid density')
    call water%calc_saturation_properties(temperatures,cp_vap=vapor,cp_liq=liquid,latent_heat=values)
    call check(vapor,cp_vap,'bulk vapor cp')
    call check(liquid,cp_liq,'bulk liquid cp')
    call check(values,latent,'bulk latent heat')
    call water%calc_saturation_properties(temperatures,drho_dT_vap=vapor,drho_dT_liq=liquid)
    call check(vapor,rho_T_vap,'bulk vapor density temperature partial')
    call check(liquid,rho_T_liq,'bulk liquid density temperature partial')
    call water%calc_saturation_properties(temperatures,drho_dP_vap=vapor,drho_dP_liq=liquid)
    call check(vapor,rho_p_vap,'bulk vapor density pressure partial')
    call check(liquid,rho_p_liq,'bulk liquid density pressure partial')
    print *, 'Independent saturation references passed'
contains
    subroutine check(actual,expected,label)
        real(real64),intent(in) :: actual(n),expected(n)
        character(*),intent(in) :: label
        integer :: j
        do j=1,n
            if (.not.ieee_is_finite(actual(j))) then
                print *, label, temperatures(j), actual(j)
                error stop 'nonfinite saturation reference result'
            end if
            ! Near critical cancellation in independent EOS evaluations requires
            ! a little more tolerance than the smooth Gibbs states.
            if (abs(actual(j)-expected(j))>2d-6*abs(expected(j))) then
                print *,label,temperatures(j),actual(j),expected(j)
                error stop 'independent saturation reference mismatch'
            end if
        end do
    end subroutine
end program
