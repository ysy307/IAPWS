# ------------------------------------------------------------------------------
# Define internal components as Object Libraries
# This avoids installing many small .a files and bundles them into the main library.
# ------------------------------------------------------------------------------

# 1. Utils: 計算補助など
add_library(iapws_utils OBJECT
    utils/kahan.F90
)

# 2. Bases: 基本型や共通定義
add_library(iapws_bases OBJECT
    bases/iapws_interface.F90
    bases/types.F90
    bases/property_helmholtz.F90
    bases/property_gibbs.F90
)

# 3. IAPWS-95: 一般水・蒸気
add_library(iapws_95 OBJECT
    iapws95/iapws95_constants.F90
    iapws95/iapws95_base.F90
    iapws95/iapws95_interface.F90
)
# 依存関係
target_link_libraries(iapws_95 PUBLIC
    iapws_utils
    iapws_bases
)

# 4. IAPWS-97: 産業用
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
# 依存関係
target_link_libraries(iapws_97  PUBLIC
    iapws_utils
    iapws_bases
)
# --- 5. IAPWS-06: 氷 (Ice Ih) ---
add_library(iapws_06 OBJECT
    iapws06/iapws06_iceIh_interface.F90
    iapws06/iapws06_iceIh_base.F90
)
target_link_libraries(iapws_06 PUBLIC  PUBLIC
    iapws_utils
    iapws_bases
)
# ------------------------------------------------------------------------------
# Main Library Target (STATIC)
# ------------------------------------------------------------------------------
add_library(iapws STATIC
    # src直下のメインファイル群
    iapws.F90
    # 分割したオブジェクトファイルを取り込む
    $<TARGET_OBJECTS:iapws_utils>
    $<TARGET_OBJECTS:iapws_bases>
    $<TARGET_OBJECTS:iapws_95>
    $<TARGET_OBJECTS:iapws_97>
    $<TARGET_OBJECTS:iapws_06>
)

# ------------------------------------------------------------------------------
# Link Libraries Configuration
# ------------------------------------------------------------------------------
target_link_libraries(iapws PRIVATE
    $<BUILD_INTERFACE:iapws_utils>
    $<BUILD_INTERFACE:iapws_bases>
    $<BUILD_INTERFACE:iapws_95>
    $<BUILD_INTERFACE:iapws_97>
    $<BUILD_INTERFACE:iapws_06>
)

# ------------------------------------------------------------------------------
# Include Directories Configuration
# ------------------------------------------------------------------------------
target_include_directories(iapws PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/iapws>
)