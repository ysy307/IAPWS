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

## Exact CPU kernel optimizations

Existing IF97 APIs also benefit without caller changes. Region-3 density solves
compute pressure and its density derivative together, using only the required
Helmholtz sums. NVHPC shares one original full coefficient evaluation to retain
bitwise compatibility. The bracketing, Newton safeguards, stopping criteria,
and EOS coefficients are unchanged. Saturation cp derivatives share Gibbs
traversals in regions 1 and 2; region 3 reuses already evaluated coefficients.
There are no lookup approximations, shared mutable caches, or MPI dependencies.

The same-API before/after benchmark is independent of any downstream application:

```sh
cmake --build build/gnu-release --target kernel_benchmark
# Build this same driver against the earlier library as well, then compare:
python3 benchmarks/compare_kernels.py --before /path/to/before/kernel_benchmark \
    --after build/gnu-release/bin/kernel_benchmark --repeats 7 --json results.json
```

It rejects differing per-sample output bits before timing, alternates build
execution order, and records medians/ranges. The optional CMake/OpenMP check
compares serial results with four workers sharing a model; its source lives
outside fpm's automatic test directory. See [the HPC kernel results](docs/hpc-kernel-results.md)
for measured gains, compiler-specific behavior, and reproduction details.

# API Document
Full API documentation is available here:
[IAPWS FORD Document](https://ysy307.github.io/IAPWS/)

## Batched IF97 evaluation

Use `calc_saturation_properties` when several quantities are needed at the same
saturation temperature. It evaluates the phase pair once, sharing saturation
pressure and Gibbs coefficients (or region-3 density solves). Requested cp
partials share the coefficient traversal in regions 1 and 2; density and
latent-heat derivatives reuse second-order coefficients. Unrequested phase cp
partials do not trigger third-order work. NVHPC retains separate Gibbs/cp
traversals to preserve the ordinary evaluator's vector-reduction rounding.
All outputs are optional; the procedure is pure elemental and accepts arrays.
Existing scalar APIs remain available with the same results.

```fortran
call water%calc_saturation_properties(T, rho_vap=rho_v, &
    density_dT_vap=drho_v_dT, cp_vap=cp_v, dcp_dT_vap=dcp_v_dT, &
    latent_heat=Lv, dlatent_heat_dT=dLv_dT)
call water%calc_cp_derivatives(T, p, dcp_dT, dcp_dp)
```

`density_dT_vap/liq` and `dcp_dT_vap/liq` are **total derivatives along
saturation**. `drho_dT_vap/liq` and `drho_dP_vap/liq` are fixed-pressure and
fixed-temperature **partials**. Additional outputs are `rho_liq`, `cp_liq`, and
the liquid counterparts of the derivative outputs above. Units, validity limits,
critical-point values, and NaN behavior match the individual methods. No saved
state or shared mutable cache is introduced.

`calc_cp_derivatives(T_in, p_in, dcp_dT, dcp_dp)` returns both ordinary partials
with one region dispatch and derivative evaluation (including one density solve
in region 3), rather than repeating the work through two individual methods.

The manual benchmark in `benchmarks/bulk_benchmark.F90` compares the six-output
saturation call and ordinary cp partials in regions 1, 2, 3, and 5. It checks each
sampled output, warms both paths, alternates execution order over seven repeats,
and reports median/min/max CPU times and matching checksums. It lives outside
`test/` so fpm auto-test discovery does not include timing runs.

```sh
cmake --build build/gnu-release --target bulk_benchmark
build/gnu-release/bin/bulk_benchmark
```

Timing is not a CTest pass criterion. These measurements compare the scalar and
bulk APIs of the same build; they do not measure an application or MPI run.
Application speedup requires adopting the new calls at consuming call sites and
measuring the full application again. For review findings, verification scope,
and measured results, see [the bulk API review](docs/bulk-api-review.md).

## State-array APIs

`calc_properties_batch` accepts state arrays and writes selected properties into
caller-owned arrays. IF97 regions 1, 2 and 5 use state-lane SIMD; IAPWS95 and ice
also expose the API through conservative evaluators. IF97 additionally provides
`calc_saturation_properties_batch`. Existing scalar APIs remain available.

```fortran
call water%calc_properties_batch(temperatures, pressures, rho=densities, cp=heat_capacities)
call water%calc_saturation_properties_batch(temperatures, latent_heat=latent_heats)
```

See the [API contract](docs/batch-api.md) and [measured performance and validation](docs/batch-performance.md).
Single-state calls should continue using the scalar API; speedups depend on the model and compiler.

# License
This project is licensed under the MIT License - see the LICENSE file for details.
