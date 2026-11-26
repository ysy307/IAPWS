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
    # The public derivative API returns quiet NaN for invalid states. Fast-math
    # assumptions can erase IEEE predicates and invalidate that contract, even
    # in callers such as the validation executables.
    add_compile_options(-ffree-line-length-none -fno-fast-math)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-fp-model=precise $<$<PLATFORM_ID:Windows>:/recursive> $<$<NOT:$<PLATFORM_ID:Windows>>:-recursive>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    add_compile_options(-Kieee)
endif()

# MSVC: specific definitions
if(MSVC)
    add_compile_definitions(WINDOWS_ENV)
endif()

# ---------------------------------------------------------
# Targets & Aliases
# ---------------------------------------------------------
set(IAPWS_BATCH_TILE_SIZE 32 CACHE STRING "States per cache-sized batch tile (4..128)")
if(NOT IAPWS_BATCH_TILE_SIZE MATCHES "^[1-9][0-9]*$" OR
   IAPWS_BATCH_TILE_SIZE LESS 4 OR IAPWS_BATCH_TILE_SIZE GREATER 128)
    message(FATAL_ERROR "IAPWS_BATCH_TILE_SIZE must be an integer between 4 and 128")
endif()

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