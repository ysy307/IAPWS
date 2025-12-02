# ------------------------------------------------------------------------------
# Internal Object Libraries
# ------------------------------------------------------------------------------

# 1. Utils
add_library(iapws_utils OBJECT
    utils/kahan.F90
)

# 2. Bases
add_library(iapws_bases OBJECT
    bases/iapws_interface.F90
    bases/types.F90
    bases/property_helmholtz.F90
    bases/property_gibbs.F90
)
target_link_libraries(iapws_bases PUBLIC iapws_utils)

# 3. IAPWS-95 (General Water/Steam)
add_library(iapws_95 OBJECT
    iapws95/iapws95_constants.F90
    iapws95/iapws95_base.F90
    iapws95/iapws95_interface.F90
)
target_link_libraries(iapws_95 PUBLIC iapws_utils iapws_bases)

# 4. IAPWS-97 (Industrial)
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
target_link_libraries(iapws_97 PUBLIC iapws_utils iapws_bases)

# 5. IAPWS-06 (Ice)
add_library(iapws_06 OBJECT
    iapws06/iapws06_iceIh_interface.F90
    iapws06/iapws06_iceIh_base.F90
)
target_link_libraries(iapws_06 PUBLIC iapws_utils iapws_bases)

# ------------------------------------------------------------------------------
# Main Library Target
# ------------------------------------------------------------------------------
add_library(iapws STATIC
    iapws.F90
    $<TARGET_OBJECTS:iapws_utils>
    $<TARGET_OBJECTS:iapws_bases>
    $<TARGET_OBJECTS:iapws_95>
    $<TARGET_OBJECTS:iapws_97>
    $<TARGET_OBJECTS:iapws_06>
)

# Set namespace for export
set_target_properties(iapws PROPERTIES EXPORT_NAME IAPWS)

# Include directories
target_include_directories(iapws PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/iapws>
)

# Link dependencies (Build Order Only)
# CRITICAL FIX: Use $<BUILD_INTERFACE:...> to hide these internal targets from the export set.
target_link_libraries(iapws PRIVATE
    $<BUILD_INTERFACE:iapws_utils>
    $<BUILD_INTERFACE:iapws_bases>
    $<BUILD_INTERFACE:iapws_95>
    $<BUILD_INTERFACE:iapws_97>
    $<BUILD_INTERFACE:iapws_06>
)

# ------------------------------------------------------------------------------
# Installation Rules
# ------------------------------------------------------------------------------
install(TARGETS iapws
    EXPORT IAPWSTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include/iapws
)

# Install Fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/iapws
    FILES_MATCHING 
    PATTERN "*.mod"
    PATTERN "*.smod"
)