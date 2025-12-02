cmake_minimum_required(VERSION 3.25)
project(IAPWS LANGUAGES Fortran VERSION 1.0.0)

enable_testing()

# ---------------------------------------------------------
# Build Output Configuration
# ---------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/.mod)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------
# Compiler Specific Settings
# ---------------------------------------------------------
# GCC: Allow long lines (essential for IAPWS formulas)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(-ffree-line-length-none)
endif()

# MSVC: specific definitions
if(MSVC)
    add_compile_definitions(WINDOWS_ENV)
endif()

# ---------------------------------------------------------
# Targets & Aliases
# ---------------------------------------------------------
add_subdirectory(src)

# Create alias for internal usage (e.g. tests or apps)
add_library(IAPWS::IAPWS ALIAS iapws)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test")
    add_subdirectory(test)
endif()

# ---------------------------------------------------------
# Package Configuration (for find_package)
# ---------------------------------------------------------
include(CMakePackageConfigHelpers)
set(CMAKE_CONFIG_INSTALL_DIR "lib/cmake/iapws")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/IAPWSTargets.cmake\")\n"
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfigVersion.cmake"
    DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
)

install(EXPORT IAPWSTargets
    FILE IAPWSTargets.cmake
    NAMESPACE IAPWS::
    DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
)