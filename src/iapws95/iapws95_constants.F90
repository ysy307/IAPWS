module iapws95_constants
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    !> Critical temperature [K]
    real(real64), parameter, public :: critical_temperature = 647.096d0
    !> Critical density [kg/m^3]
    real(real64), parameter, public :: critical_density = 322.0d0

    !> Specific gas constant for water [J/(kg·K)]
    real(real64), parameter, public :: specific_gas_constant_water = 0.46151805d3

    !=================================================================================
    ! IAPWS-95 ideal helmholtz energy equation coefficients
    !=================================================================================
    real(real64), parameter, public :: n0_log(2) = [1.0d0, 3.00632d0]
    integer(int32), parameter, public :: pow(2) = [0, 1]
    real(real64), parameter, public :: n0_pow(2) = [-8.3204464837497d0, 6.6832105275932d0]

    real(real64), parameter, public :: n0_exp(5) = [ &
                                       0.012436d0, 0.97315d0, 1.27950d0, &
                                       0.96956d0, 0.24873d0]
    real(real64), parameter, public :: g0(5) = [ &
                                       1.28728967d0, 3.53734222d0, 7.74073708d0, &
                                       9.24437796d0, 27.5075105d0]

    !=================================================================================
    ! End of IAPWS-95 residual helmholtz energy equation coefficients
    !=================================================================================
    real(real64), parameter, public :: nr1(7) = [ &
                                       0.12533547935523d-1, 0.78957634722828d1, -0.87803203303561d1, &
                                       0.31802509345418d0, -0.26145533859358d0, -0.78199751687981d-2, &
                                       0.88089493102134d-2]
    integer(int32), parameter, public :: d1(7) = [1, 1, 1, 2, 2, 3, 4]
    real(real64), parameter, public :: t1(7) = [-0.5d0, 0.875d0, 1.0d0, 0.5d0, 0.75d0, 0.375d0, 1.0d0]

    real(real64), parameter, public :: nr2(44) = [ &
                                       -0.66856572307965d0, 0.20433810950965d0, -0.66212605039687d-4, &
                                       -0.19232721156002d0, -0.25709043003438d0, 0.16074868486251d0, &
                                       -0.40092828925807d-1, 0.39343422603254d-6, -0.75941377088144d-5, &
                                       0.56250979351888d-3, -0.15608652257135d-4, 0.11537996422951d-8, &
                                       0.36582165144204d-6, -0.13251180074668d-11, -0.62639586912454d-9, &
                                       -0.10793600908932d0, 0.17611491008752d-1, 0.22132295167546d0, &
                                       -0.40247669763528d0, 0.58083399985759d0, 0.49969146990806d-2, &
                                       -0.31358700712549d-1, -0.74315929710341d0, 0.47807329915480d0, &
                                       0.20527940895948d-1, -0.13636435110343d0, 0.14180634400617d-1, &
                                       0.83326504880713d-2, -0.29052336009585d-1, 0.38615085574206d-1, &
                                       -0.20393486513704d-1, -0.16554050063734d-2, 0.19955571979541d-2, &
                                       0.15870308324157d-3, -0.16388568342530d-4, 0.43613615723811d-1, &
                                       0.34994005463765d-1, -0.76788197844621d-1, 0.22446277332006d-1, &
                                       -0.62689710414685d-4, -0.55711118565645d-9, -0.19905718354408d0, &
                                       0.31777497330738d0, -0.11841182425981d0]
    integer(int32), parameter, public :: c2(44) = [ &
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                         1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, &
                                         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                                         2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]
    integer(int32), parameter, public :: d2(44) = [ &
                                         1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, &
                                         10, 11, 13, 15, 1, 2, 2, 2, 3, 4, &
                                         4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, &
                                         10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
    integer(int32), parameter, public :: t2(44) = [ &
                                         4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, &
                                         7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, &
                                         16, 22, 23, 23, 10, &
                                         50, 44, 46, 50]
    integer(int32), parameter, public :: gamma2(44) = [ &
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    real(real64), parameter, public :: nr3(3) = [ &
                                       -0.31306260323435d2, 0.31546140237781d2, -0.25213154341695d4]
    integer(int32), parameter, public :: d3(3) = [3, 3, 3]
    integer(int32), parameter, public :: t3(3) = [0, 1, 4]
    real(real64), parameter, public :: alfa3(3) = [ &
                                       20.0d0, 20.0d0, 20.0d0]
    real(real64), parameter, public :: beta3(3) = [ &
                                       150.0d0, 150.0d0, 250.0d0]
    real(real64), parameter, public :: gamma3(3) = [ &
                                       1.21d0, 1.21d0, 1.25d0]
    real(real64), parameter, public :: epsilon3(3) = [ &
                                       1.0d0, 1.0d0, 1.0d0]

    real(real64), parameter, public :: nr4(2) = [ &
                                       -0.14874640856724d0, 0.31806110878444d0]
    real(real64), parameter, public :: a4(2) = [3.5d0, 3.5d0]
    real(real64), parameter, public :: b4(2) = [0.85d0, 0.95d0]
    real(real64), parameter, public :: B(2) = [0.2d0, 0.2d0]
    integer(int32), parameter, public :: C(2) = [28, 32]
    integer(int32), parameter, public :: D(2) = [700, 800]
    real(real64), parameter, public :: A(2) = [0.32d0, 0.32d0]
    real(real64), parameter, public :: beta4(2) = [0.3d0, 0.3d0]
end module iapws95_constants
