# IAPWS: Fortran Library for Water and Steam Properties

This library provides a modern Fortran implementation of the thermodynamic properties of water and steam, based on the formulations released by the International Association for the Properties of Water and Steam (IAPWS).

It is designed for scientific simulations and industrial applications requiring high-precision calculation of fluid properties.

## Features
The library implements the following IAPWS standards:
- **IAPWS-95**: The formulation for general and scientific use (high accuracy).
- **IAPWS-IF97**: The industrial formulation for calculating the properties of water and steam (optimized for computational speed).
- **IAPWS-06**: The formulation for the properties of Ice Ih. 

## Technical Highlights
- Written in Modern Fortran (F2008+ standards).
- Supports multiple build systems: CMake (Recommended) and Fortran Package Manager (fpm).
- Cross-platform compatibility (Linux, Windows, macOS).

## Prerequisites
To build this library, you need a Fortran compiler and a build system.
- Fortran Compiler:
    - GNU Fortran (**gfortran**)
    - Intel OneAPI Compiler (**ifx**)
    - NVIDIA HPC SDK (**nvfortran**)
- Build Tools (either one):
    - CMake (Version 3.25 or later) - Recommended
    - Fortran Package Manager (fpm)
    
## Build and Installation
CMake is the primary and recommended build system for installing and integrating this library.
### Method 1:
Using CMake (Recommended)This method is best for system-wide installation or integration into C++/mixed projects.

#### Basic Build:
```sh
cmake -S . -B build
cmake --build build
```
Using Presets (Cross-Platform):This project provides CMake presets for different compilers (`gcc`, `intel`, `nvidia`, `windows-msvc`).
```sh
# List available presets
cmake --list-presets

# Configure (e.g., using GCC Release settings)
cmake --preset gcc-release

# Build
cmake --build --preset build-gcc-release

# Test
ctest --preset test-gcc-release
```
#### Installation:
```sh
cmake --install build --prefix /usr/local
```

### Method 2: Using the Helper Script (For Developers)
A helper script FPMBuild.sh is provided to simplify the build and test process using fpm with predefined compiler flags. This is useful for rapid development and testing.

#### Usage:
```sh
./FPMBuild.sh [compiler] [build_type]
```
Arguments:
- `compiler`: `gcc`, `intel`, `nvidia`, or `clean`
- `build_type`: `debug` or `release`

#### Examples:
```sh
# Build and test using gfortran in debug mode
./FPMBuild.sh gcc debug

# Build and test using Intel ifx in release mode
./FPMBuild.sh intel release

# Clean build artifacts
./FPMBuild.sh clean
```

### Method 3: Using Fortran Package Manager (fpm)
If you prefer using standard fpm commands for development:
```sh
# Build the project
fpm build

# Run tests
fpm test

# Build with release optimizations
fpm build --profile release
```

## Usage in Your Project
### Using CMake
In your CMakeLists.txt, use `find_package`:
```make
find_package(IAPWS REQUIRED)
add_executable(my_app main.f90)
target_link_libraries(my_app PRIVATE IAPWS::IAPWS)
```

### Using fpm
Add the following dependency to your fpm.toml:
```toml
[dependencies]
IAPWS = { git = "https://github.com/ysy307/IAPWS.git" }
```
