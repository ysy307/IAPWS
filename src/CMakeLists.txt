# ライブラリの定義
add_library(iapws STATIC
    iapws95/iapws95_interface.F90
)

# モジュールファイルの参照パスを公開
target_include_directories(iapws PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)