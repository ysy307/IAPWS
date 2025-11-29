# ライブラリの定義
add_library(iapws STATIC
    iapws_interface.F90
    bases/properties.F90
    iapws95/iapws95_interface.F90
    iapws95/iapws95_constants.F90
    iapws95/iapws95_base.F90
    iapws.F90
)

# モジュールファイルの参照パスを公開
target_include_directories(iapws PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)