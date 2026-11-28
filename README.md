# IAPWS: Fortran Library for Water and Steam Properties

[![CI Status](https://github.com/ysy307/IAPWS/actions/workflows/verification.yml/badge.svg)](https://github.com/ysy307/IAPWS/actions)
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
./scripts/FPMBuild.sh [compiler] [build_type]
```
Arguments:
- `compiler`: `gcc`, `intel`, `nvidia`, or `clean`
- `build_type`: `debug` or `release`

#### Examples:
```sh
# Build and test using gfortran in debug mode
./scripts/FPMBuild.sh gcc debug

# Build and test using Intel ifx in release mode
./scripts/FPMBuild.sh intel release

# Clean build artifacts
./scripts/FPMBuild.sh clean
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

## Fortran static analysis

IAPWS includes the FortranCheck tool developed in FTCMS. Build it with Rust 1.88 or later, Cargo, and a C compiler. The CMake integration supports single-configuration Ninja or Makefiles with `gfortran` or `ifx`. [FortranCheck documentation](tools/fortran-check/README.md)

```sh
cargo install --path tools/fortran-check --locked
cmake --preset gnu-debug -B build/gnu-check -DIAPWS_ENABLE_FORTRAN_CHECK=ON
cmake --build build/gnu-check --target fortran_check
```

Set `FORTRAN_CHECK_EXECUTABLE` at configure time to select an executable outside `PATH`. In this build directory, CMake runs the check before compiling and uses the flags and definitions from its compilation database to analyze `src/`.

The policy permits `INTENT(OUT)` for result arguments and preserves published data fields and procedure names. Explicit complexity ceilings cover existing numerical kernels and public model types. Dependency errors and incomplete analysis fail the check. [IAPWS policy](fortran-check.toml)

The reports in `build/gnu-check/fortran-check/` include `diagnostics.sarif`, `metrics.json`, and `dependencies.dot`. The GitHub Actions verification workflow uploads these reports. Documentation publishes from `main` after the static check, GNU Debug and Release tests, and documentation build succeed.

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

## Unified pure-water API

Construct a `type_water` with `new_water(fluid=IF97)` for industrial fluid
properties or `new_water(fluid=IAPWS95)` for scientific fluid and ice equilibrium
calculations. Inputs and outputs use SI units.

```fortran
use iso_fortran_env, only: real64
use iapws, only: type_water, type_water_state, type_water_saturation, &
                 new_water, IF97, IAPWS95, LIQUID, WATER_OK, WATER_METASTABLE

type(type_water) :: water, equilibrium
type(type_water_state) :: state
type(type_water_saturation) :: sat
real(real64) :: mu, k

water = new_water(fluid=IF97)
call water%state_tp(400.0_real64, 1.0e5_real64, state)
if (state%status == WATER_OK) call water%transport(state, mu=mu, conductivity=k)
call water%saturation_p(1.0e5_real64, sat)

equilibrium = new_water(fluid=IAPWS95)
call equilibrium%state_ph(1.0e5_real64, -150000.0_real64, state)
! For a successful ice-liquid mixture, state%solid_fraction is its ice mass fraction.
call equilibrium%state_tp(250.0_real64, 1.0e5_real64, state, phase=LIQUID, allow_metastable=.true.)
! G12-15 return values have WATER_METASTABLE status.
```

- `state_tp`, `state_ph`, `state_ps`, `saturation_t` and `saturation_p` are
  pure elemental procedures and accept conformable arrays. A state reports its
  `phase`, actual `formulation`, `status` and `properties`. Coexisting phases are
  marked in `coexistence`, with bit positions given by the phase constants.
  A `T,p` coexistence query can select one endpoint with `phase=`.
- R14-08 supplies the melting and sublimation boundaries. `state_tp` returns
  Ice Ih from IAPWS-06 even when the selected fluid equation is IF97. Higher
  pressure ice can be detected, with its thermodynamic properties reported as
  unsupported. `melting_pressure(T, solid)` also provides the individual R14-08
  melting curves through ice VII.
- Solid-phase `p,h` and `p,s` evaluations use the IAPWS-95 constructor. They
  combine IAPWS-95 (or G9-12 for low-temperature vapor) with IAPWS-06 at the
  R14-08 boundary. Phase fractions at the triple point require additional
  information. IF97's optional `extend_cold_liquid=.true.` applies to
  `state_tp` and uses IAPWS-95 below the IF97 lower temperature limit in
  stable liquid states.
- Explicit phase and metastable requests enable G12-15 supercooled liquid
  and the IF97 equation for metastable steam through 10 MPa. They report
  failure outside each equation's recommended range. G9-12 extends scientific
  vapor calculations down to 50 K in the stable vapor region.
- Wet IF97 or IAPWS-95 states return vapor mass fraction in `quality`. Mixture
  density follows from the mass-weighted specific volume. Single-phase
  response functions such as `cp` and sound speed are NaN for mixtures.
- `transport` returns dynamic viscosity in Pa·s and thermal conductivity
  in W/(m·K) for supported fluid states. Saturation results include
  `sigma` in N/m. When using the standalone transport functions with
  only temperature and density, the result contains the dilute-gas
  and finite-density terms.
- Near a join where two IF97 regions produce the same enthalpy or entropy
  at a fixed pressure, the inverse API can return `WATER_AMBIGUOUS`. Pass
  `region=` to select the required IF97 branch.
- At the exact critical point the common API returns `WATER_SINGULAR`
  with critical coordinates and density. A failed saturation solve very
  close to the critical point returns `WATER_NO_CONVERGENCE`.

### CFD and thermal-hydraulic solver interface

```fortran
use iso_fortran_env, only: real64
use iapws, only: type_water, type_water_state, type_water_jacobian, &
                 new_water, IAPWS95, WATER_OK, JAC_RHO, JAC_T

type(type_water) :: water
type(type_water_state) :: state
type(type_water_jacobian) :: jac
real(real64) :: p, h, previous_T, drho_dp, drho_dh, dT_dh

water = new_water(fluid=IAPWS95)
p = 3.0e6_real64
h = 5.3e5_real64
previous_T = 400.0_real64
call water%state_ph_checked(p, h, state, guess_T=previous_T, jac=jac)
if (state%status == WATER_OK .and. jac%ph_available) then
    drho_dp = jac%ph(JAC_RHO, 1)
    drho_dh = jac%ph(JAC_RHO, 2)
    dT_dh = jac%ph(JAC_T, 2)
end if
```

- `state_tp_checked`, `state_ph_checked`, and `state_ps_checked` return
  `WATER_OK` for successful equilibrium states and `WATER_METASTABLE` for
  successful explicit metastable queries. On failure the property record
  contains finite zero placeholders and the status indicates the reason.
  Check the status before using any property.
- `jac%tp`, `jac%ph`, and `jac%ps` have columns ordered as their respective
  input pairs, with the other input held fixed. Rows are `JAC_T`, `JAC_P`,
  `JAC_RHO`, `JAC_V`, `JAC_H`, `JAC_S`, and `JAC_U`. Check the
  corresponding availability flag. The matrices use SI units.
- For equilibrium two-phase states the `ph` and `ps` matrices differentiate
  the same boundary and mass-weighted properties used by the inverse.
  Temperature and pressure alone leave the phase fraction undetermined.
  The checked API fills undefined single-phase quantities such as `cp` and
  sound speed with zeros in these states. Use the phase and availability flags
  to interpret the outputs. At a phase endpoint `phase_endpoint` indicates
  that the derivative follows the returned branch. IF97 region joins retain the
  published piecewise equations and the returned `properties%region_id`
  identifies the derivative's region.
- Pass the previous cell iterate as `guess_T=` for `ph` and `ps`, and
  `guess_rho=` for `tp` with IAPWS-95 or IF97 region 3. Hints are optional.
  A rejected hint falls back to the bracketed solver. A known IF97 region
  can also be passed with `region=` to specify the inverse branch.
- The checked APIs and Jacobian procedure are pure elemental. They can be
  called with conformable arrays or inside an OpenMP loop sharing the
  constructed `type_water`. The caller supplies the guesses independently
  for each cell. The façade owns the exact EOS evaluators privately, and the
  public state and Jacobian types carry the solver results.

The verification fixtures include IAPWS R6-95(2018) Table 7, R7-97(2012)
region check values, R10-06(2009) Table 6 states, R14-08 Table 3,
G12-15 Table 5, G9-12 Table 3, R12-08 and R15-11 Tables 4 and 5, and
SR5-05 Tables 5 and 12. The published values are rounded. The tests also compare unrounded
results with Python `iapws==1.5.5` if installed. IEEE arithmetic is required.
[Official Releases and Guidelines](https://iapws.org/technical-guidance/release)

```sh
python3 -m pip install iapws==1.5.5
cmake --workflow --preset workflow-gnu-release
cmake --build --preset build-gnu-release --target water_performance
build/gnu-release/bin/water_performance
OMP_NUM_THREADS=4 OMP_PROC_BIND=close build/gnu-release/bin/water_performance --parallel
```

The performance driver reports twelve IF97 region-3 densities and the time
for 120,000 density evaluations in seconds and ns/point. Compile the same driver against an earlier
revision with the same compiler and flags for a comparison. The accelerated
solver uses SR5-05 as a seed, refines against the forward pressure equation,
and retains the bracketed solver near the critical point and as a fallback.
The refinement groups the pressure coefficients at the fixed temperature,
then evaluates pressure and its density derivative with Horner's method.
The SR5-05 subregion fits evaluate only their required exponents. The first
refinement uses a guarded Halley step and subsequent steps use Newton updates.
The benchmark checks every timed result against its scalar density. Its
explicit parallel mode uses OpenMP across independent evaluations.

## External water comparison

The optional CoolProp test compares like formulations in SI units.
`HEOS::Water` is used for IAPWS-95 and `IF97::Water` for IF97. It checks
single-phase density, heat capacity, enthalpy, saturation pressure and
endpoint densities. The HEOS comparison also checks three analytic
`(p,h)` Jacobian entries. Each native single-phase state is inverted back
to temperature with a tolerance of `1e-8 K`. The exact critical point must
return `WATER_SINGULAR`. Ice Ih is covered by the existing IAPWS-06
reference fixtures and Python `iapws` comparison.

```sh
python3 -m pip install CoolProp==8.1.0
cmake --build --preset build-gnu-release --target edge_states
ctest --test-dir build/gnu-release -R coolprop_reference -V
# Also compare REFPROP when the native library and water data are installed.
python3 validations/compare_water.py build/gnu-release/bin/edge_states --refprop
```

Observed on GNU Release with CoolProp 8.1.0. Values are maximum relative
differences across the sampled grid. The test reports enthalpy separately
in J/kg because it approaches zero near the triple point.

| Comparison | Triple vicinity saturation density | Critical saturation density, 0.01–1 K below | Critical saturation density, 0.0001–0.001 K below |
|---|---:|---:|---:|
| IAPWS-95 / HEOS | 1.08e-14 | 8.78e-10 | 1.97e-7 |
| IF97 / CoolProp IF97 | 2.13e-15 | 3.65e-7 | 4.72e-4 |

The sampled-grid relative tolerances for ordinary, triple vicinity, critical
and inner critical cases are respectively `5e-12`, `5e-9`, `5e-8` and
`5e-6` for HEOS, and `5e-12`, `5e-12`, `5e-6` and `5e-3` for IF97.
These are regression limits for this grid. Three IF97 vapor-side `T,p`
points at 605.54 Pa fall below the CoolProp IF97 pressure limit. The test
still requires successful native evaluation and inverse at those points.
See the [CoolProp IF97 backend documentation](https://coolprop.org/fluid_properties/IF97.html).

The native serial density benchmark on an Intel Core i7-14700 with GNU Release
took **0.0146 s for 120,000 points (122 ns/point)**, the median of five runs
(range 116–126 ns/point). The grid repeats twelve fixed points. This measures
the native density call plus the benchmark's checksum and verification work.

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
bitwise compatibility. The EOS coefficients are unchanged. Saturation cp derivatives share Gibbs
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
