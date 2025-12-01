cmake_minimum_required(VERSION 3.25)
project(IAPWS LANGUAGES Fortran)
enable_testing()
# ---------------------------------------------------------
# Settings for build output directories
# ---------------------------------------------------------
# Module files -> <build>/.mod/
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/.mod)
# Static libraries (.a) -> <build>/lib/
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# Executable files -> <build>/bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------
# Reading subdirectories
# ---------------------------------------------------------
add_subdirectory(src)
add_subdirectory(app)

# ---------------------------------------------------------
# Installation settings
# ---------------------------------------------------------
install(TARGETS iapws
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
        DESTINATION include
        FILES_MATCHING PATTERN "*.mod"
)
