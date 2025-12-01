# --- 1. 内部モジュール (OBJECTライブラリ) の定義 ---
# ※ ここでのパスはすべて src/ からの相対パスになります

# Utils: 計算補助など
add_library(iapws_utils OBJECT
    utils/kahan.F90
)

# Bases: 基本型や共通定義
add_library(iapws_bases OBJECT
    iapws_interface.F90
    bases/types.F90
    bases/property_helmholtz.F90
    bases/property_gibbs.F90
)
# 依存関係 (Bases は Utils を使う)
target_link_libraries(iapws_bases PUBLIC iapws_utils)

# IAPWS-95: 一般水・蒸気
add_library(iapws_95 OBJECT
    iapws95/iapws95_constants.F90
    iapws95/iapws95_base.F90
    iapws95/iapws95_interface.F90
)
# 依存関係 (95 は Bases を使う)
target_link_libraries(iapws_95 PUBLIC iapws_bases)

# IAPWS-97: 産業用
add_library(iapws_97 OBJECT
    iapws97/iapws97_base.F90
    iapws97/iapws97_auxiliary.F90
    iapws97/iapws97_region1.F90
    iapws97/iapws97_region2.F90
    iapws97/iapws97_region3.F90
    iapws97/iapws97_region4.F90
    iapws97/iapws97_region5.F90
    iapws97/iapws97_interface.F90
)
# 依存関係 (97 は Bases を使う)
target_link_libraries(iapws_97 PUBLIC iapws_bases)


# --- 2. 最終ライブラリの定義 (STATIC) ---
# ここで全てのオブジェクトを結合して libiapws.a を作成

add_library(iapws STATIC
    # src直下のファイル
    iapws.F90
    # 分割したオブジェクトファイルを取り込む
    $<TARGET_OBJECTS:iapws_utils>
    $<TARGET_OBJECTS:iapws_bases>
    $<TARGET_OBJECTS:iapws_95>
    $<TARGET_OBJECTS:iapws_97>
)

# Fortranのモジュール参照(.mod)の順序解決のためにリンク設定を行う
target_link_libraries(iapws PRIVATE
    iapws_utils
    iapws_bases
    iapws_95
    iapws_97
)

# --- 3. インターフェースの公開 ---
# このライブラリを使う側 (appなど) に .mod ファイルの場所を教える
target_include_directories(iapws PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)