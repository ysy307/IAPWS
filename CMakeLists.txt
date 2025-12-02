cmake_minimum_required(VERSION 3.25)
project(IAPWS LANGUAGES Fortran VERSION 1.0.0)
enable_testing()

# ---------------------------------------------------------
# Configuration for Build Output Directories
# ---------------------------------------------------------
# Specify the directory where Fortran module files (.mod) are placed during the build.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/.mod)
# Specify the directory where static libraries (.a) are placed during the build.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# Specify the directory where executable files are placed during the build.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------
# Process Subdirectories
# ---------------------------------------------------------
add_subdirectory(src)
add_subdirectory(app)

# ---------------------------------------------------------
# Installation Configuration
# ---------------------------------------------------------

# 1. Install Library Targets and Export Information
# Install the static library to <install_directory>/lib/
# Install the executables to <install_directory>/bin/
# Register the include directory as <install_directory>/include/iapws
install(TARGETS iapws
    EXPORT IAPWSTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include/iapws
)

# 2. Install Fortran Module Files
# Copy all .mod and .smod files to <install_directory>/include/iapws/
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/iapws
    FILES_MATCHING 
    PATTERN "*.mod"
    PATTERN "*.smod"
)

# ---------------------------------------------------------
# CMake Package Configuration File Generation
# ---------------------------------------------------------
include(CMakePackageConfigHelpers)

# Define the installation directory for CMake configuration files.
# Path: <install_directory>/lib/cmake/iapws/
set(CMAKE_CONFIG_INSTALL_DIR "lib/cmake/iapws")

# Generate the Version Configuration File (IAPWSConfigVersion.cmake)
# This file allows find_package() to check for version compatibility.
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Generate the Main Configuration File (IAPWSConfig.cmake)
# This file includes the exported targets file so that downstream projects can import the targets.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/IAPWSTargets.cmake\")\n"
)

# Install the Generated Configuration Files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/IAPWSConfigVersion.cmake"
    DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
)

# Install the Exported Targets File (IAPWSTargets.cmake)
# This defines the imported targets (e.g., IAPWS::iapws) for other projects to use.
install(EXPORT IAPWSTargets
    FILE IAPWSTargets.cmake
    NAMESPACE IAPWS::
    DESTINATION ${CMAKE_CONFIG_INSTALL_DIR}
)