# ライブラリの定義
add_library(iapws STATIC
    iapws_interface.F90
    bases/property_helmholtz.F90
    bases/property_gibbs.F90
    bases/types.F90
    utils/kahan.F90
    iapws95/iapws95_interface.F90
    iapws95/iapws95_constants.F90
    iapws95/iapws95_base.F90
    iapws97/iapws97_interface.F90
    iapws97/iapws97_auxiliary.F90
    iapws97/iapws97_region1.F90
    iapws97/iapws97_region2.F90
    iapws97/iapws97_region3.F90
    iapws97/iapws97_region4.F90
    iapws97/iapws97_region5.F90
    iapws.F90
)

# モジュールファイルの参照パスを公開
target_include_directories(iapws PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)