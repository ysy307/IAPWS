cmake_minimum_required(VERSION 3.25)
project(IAPWS LANGUAGES Fortran)
enable_testing()
# ---------------------------------------------------------
# Settings for build output directories
# ---------------------------------------------------------
# Module files -> <build>/.mod/
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/.mod)
# Static libraries (.a) -> lib/
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
# Executable files -> bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# ---------------------------------------------------------
# Reading subdirectories
# ---------------------------------------------------------
add_subdirectory(src)
add_subdirectory(app)