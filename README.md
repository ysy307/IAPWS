# IAPWS: Fortran Library for Water and Steam Properties

[![CI Status](https://github.com/ysy307/IAPWS/actions/workflows/documentation.yml/badge.svg)](https://github.com/ysy307/IAPWS/actions)
[![License](https://img.shields.io/github/license/ysy307/IAPWS)](LICENSE)

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
- Cross-platform compatibility (Linux, Windows).

## Prerequisites
To build this library, you need a Fortran compiler and a build system.
- Fortran Compiler:
    - GNU Fortran (`gfortran`)
    - Intel OneAPI Compiler (`ifx`)
    - NVIDIA HPC SDK (`nvfortran`)
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
Using Presets (Cross-Platform):This project provides CMake presets for different compilers (`gnu`, `intel`, `nvidia`, `windows-msvc`).
```sh
# List available presets
cmake --list-presets

# Configure (e.g., using GCC Release settings)
cmake --preset gnu-release

# Build
cmake --build --preset build-gnu-release

# Test
ctest --preset test-gnu-release
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

## Analytic derivatives

Initialize the model before calling these pure elemental `real64` methods. Scalar
and conformable array arguments are supported. Temperatures are in K and
pressures in Pa.

| Model | Method | Derivative / output units |
| --- | --- | --- |
| IF97 and Ice Ih | `calc_dcp_dT(T_in, p_in, dcp_dT)` | At fixed pressure, J/(kg K²) |
| IF97 and Ice Ih | `calc_dcp_dp(T_in, p_in, dcp_dp)` | At fixed temperature, J/(kg K Pa) |
| IF97 | `calc_saturation_dcp_dT(T_in, dcp_dT_vap, dcp_dT_liq)` | Along saturation, J/(kg K²) |
| IF97 | `calc_latent_heat_dT(T_in, dlatent_heat_dT)` | Along saturation, J/(kg K) |
| IF97 | `calc_saturation_density_dT(T_in, drho_dT_vap, drho_dT_liq)` | Along saturation, kg/(m³ K) |

Fortran names are case insensitive; the Ice pressure argument is named `P_in`.
Both phase outputs are optional and ordered vapor, liquid. Saturation derivatives
include the analytic pressure slope `dp_sat/dT`. The existing
`calc_saturation_drho_dT` and `calc_saturation_drho_dP` methods retain their
**partial** derivative meanings at fixed pressure and temperature, respectively.

IF97 ordinary derivatives follow the existing region dispatcher: 243.15–2273.15 K,
0 < p ≤ 100 MPa, with p ≤ 50 MPa above 1073.15 K. The sub-273.15 K range preserves
the library's existing extrapolation convention. Ice heat-capacity derivatives
evaluate the same Gibbs expression as `calc_cp` for finite T and p, including
extrapolation outside 0 < T ≤ 273.16 K and 0 ≤ p ≤ 210 MPa (also at T ≤ 0
and p < 0). Extrapolation does not imply physical validity or guaranteed accuracy;
nonfinite inputs return quiet NaN. Saturation derivatives accept
273.15 ≤ T < 647.096 K, selecting regions 1/2 through 623.15 K and distinct stable
region-3 vapor/liquid roots above it. At region joins the derivative belongs to
the selected branch; the piecewise fits need not have matching derivatives.

Invalid/nonfinite inputs and critical-point singularities return quiet NaN for
the affected element. Extremely close to the critical point, the independently
fitted saturation pressure and region-3 equation of state can cease to bracket
both stable phases; saturation methods then return NaN for both phases rather
than an unstable root. At the critical point, saturation densities are 322 kg/m³
and latent heat is zero; saturation heat capacity and derivatives return NaN.
CMake builds preserve IEEE behavior even with the supplied Release presets;
other build systems must also avoid fast-math assumptions that erase NaN checks.

```fortran
call water%initialize()
call water%calc_dcp_dT(500.0_real64, 3.0e6_real64, cp_temperature_slope)
call water%calc_saturation_density_dT(640.0_real64, &
    drho_dT_vap=vapor_density_slope, drho_dT_liq=liquid_density_slope)
```

The CTest suite checks analytic results against numerical differentiation at
multiple steps, phase references, region joins, optional outputs, and invalid
array elements. The reference-data validator fails on mismatches or nonfinite
values and writes logs under the build directory. Run it directly with
`build/gnu-debug/bin/iapws_test /absolute/path/to/validations` from a scratch
working directory. Python 3, when available, enables an additional validator
regression test; the library has no Python dependency.

# API Document
Full API documentation is available here:
[IAPWS FORD Document](https://ysy307.github.io/IAPWS/)

# License
This project is licensed under the MIT License - see the LICENSE file for details.
