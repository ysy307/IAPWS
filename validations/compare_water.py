#!/usr/bin/env python3
"""Compare matching water formulations using CoolProp's low-level SI API."""

import argparse
from collections import defaultdict
import math
import subprocess
import sys

try:
    import CoolProp.CoolProp as cp
except ModuleNotFoundError as error:
    if error.name != "CoolProp":
        raise
    print("Install CoolProp to run the external water comparison.", file=sys.stderr)
    sys.exit(77)

TRIPLE_T = 273.16
TRIPLE_P = 611.657
CRITICAL_T = 647.096
CRITICAL_P = 22.064e6

# Limits apply to the sampled grid. Enthalpy is absolute because h is near zero
# at the triple point. The other property and Jacobian errors are relative.
TOLERANCES = {
    "HEOS": {
        "ordinary": (5e-12, 1e-5),
        "triple": (5e-9, 1e-5),
        "critical": (5e-8, 1e-3),
        "critical-inner": (5e-6, 1e-2),
    },
    "IF97": {
        "ordinary": (5e-12, 1e-5),
        "triple": (5e-12, 1e-5),
        "critical": (5e-6, 0.1),
        "critical-inner": (5e-3, 5.0),
    },
}
ROUNDTRIP_TOLERANCE_K = 1e-8


def cases_for(formulation):
    for temperature, pressure in [(300, 3e6), (450, 1e6), (650, 25e6), (700, 3e6)]:
        yield formulation, "tp", "ordinary", temperature, pressure
    for offset in [1e-5, 1e-3, 0.1]:
        temperature = TRIPLE_T + offset
        for pressure_factor in [0.99, 1.01]:
            yield formulation, "tp", "triple", temperature, TRIPLE_P * pressure_factor
        yield formulation, "sat", "triple", temperature, 0
    yield formulation, "sat", "triple", TRIPLE_T, 0

    for offset in [1, 1e-1, 1e-2, 1e-3, 1e-4]:
        group = "critical" if offset >= 1e-2 else "critical-inner"
        yield formulation, "sat", group, CRITICAL_T - offset, 0
        yield formulation, "tp", group, CRITICAL_T + offset, CRITICAL_P
    yield formulation, "tp", "critical-exact", CRITICAL_T, CRITICAL_P


def reference(state, operation, temperature, pressure, derivatives):
    if operation == "sat":
        state.update(cp.QT_INPUTS, 0, temperature)
        saturation_pressure, liquid_density = state.p(), state.rhomass()
        state.update(cp.QT_INPUTS, 1, temperature)
        return dict(p=saturation_pressure, rho_l=liquid_density, rho_v=state.rhomass())

    state.update(cp.PT_INPUTS, pressure, temperature)
    result = dict(rho=state.rhomass(), cp=state.cpmass(), h=state.hmass())
    if derivatives:
        result["rho_p|h"] = state.first_partial_deriv(cp.iDmass, cp.iP, cp.iHmass)
        result["rho_h|p"] = state.first_partial_deriv(cp.iDmass, cp.iHmass, cp.iP)
        result["T_h|p"] = state.first_partial_deriv(cp.iT, cp.iHmass, cp.iP)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edge_states", help="path to the compiled Fortran comparison driver")
    parser.add_argument("--refprop", action="store_true", help="also compare an installed REFPROP backend")
    args = parser.parse_args()

    backends = [(95, "HEOS"), (97, "IF97")]
    if args.refprop:
        backends.append((95, "REFPROP"))
    print(f"CoolProp {cp.get_global_param_string('version')}")
    failures = []
    for formulation, backend in backends:
        cases = list(cases_for(formulation))
        try:
            state = cp.AbstractState(backend, "Water")
        except (ValueError, OSError) as error:
            parser.error(f"{backend}::Water could not be initialized: {error}")
        rows = subprocess.run(
            [args.edge_states],
            input="".join(f"{model} {operation} {temperature:.17g} {pressure:.17g}\n"
                          for model, operation, _, temperature, pressure in cases),
            text=True, capture_output=True, check=True,
        ).stdout.splitlines()
        maxima = defaultdict(lambda: defaultdict(float))
        counts = defaultdict(lambda: [0, 0])
        if len(rows) != len(cases):
            raise RuntimeError("unexpected Fortran row count")
        print(f"\nIAPWS-{formulation} vs {backend}::Water")
        for (_, operation, group, temperature, pressure), row in zip(cases, rows):
            counts[group][0] += 1
            columns = row.split()
            location = f"{operation} T={temperature:.9f} K p={pressure:.9g} Pa"
            status = int(columns[0])
            if group == "critical-exact":
                print(f"  {location}: Fortran status {status} (expected singular=7)")
                if status != 7:
                    raise AssertionError(f"exact critical point returned status {status}")
                continue
            if status:
                failures.append(f"{backend} {location}: Fortran status {status}")
                continue
            if operation == "tp":
                inverse_status = int(columns[7])
                if inverse_status:
                    failures.append(f"{backend} {location}: Fortran p,h inverse status {inverse_status}")
                else:
                    error = abs(float(columns[8]) - temperature)
                    maxima[group]["p,h roundtrip K"] = max(maxima[group]["p,h roundtrip K"], error)
            try:
                expected = reference(state, operation, temperature, pressure, backend == "HEOS")
            except (ValueError, IndexError) as error:
                message = f"{location}: {backend} rejected input ({error})"
                if backend != "IF97" or group != "triple" or operation != "tp" or pressure >= TRIPLE_P:
                    failures.append(message)
                else:
                    print(f"  {message}")
                continue
            counts[group][1] += 1
            for col, (field, target) in enumerate(expected.items(), start=1):
                actual = float(columns[col])
                difference = abs(actual - target)
                if not math.isfinite(difference):
                    raise AssertionError(f"{location}: nonfinite comparison for {field}")
                if field == "h":
                    maxima[group]["h absolute J/kg"] = max(maxima[group]["h absolute J/kg"], difference)
                else:
                    relative = difference / abs(target) if target else difference
                    maxima[group][field] = max(maxima[group][field], relative)
        for group, (attempted, compared) in counts.items():
            if group == "critical-exact":
                continue
            for field, error in maxima[group].items():
                if field == "p,h roundtrip K":
                    limit = ROUNDTRIP_TOLERANCE_K
                elif backend in TOLERANCES:
                    relative_limit, enthalpy_limit = TOLERANCES[backend][group]
                    limit = enthalpy_limit if field == "h absolute J/kg" else relative_limit
                else:
                    continue
                if error > limit:
                    failures.append(f"{backend} {group}: {field}={error:.3e} exceeds {limit:.3e}")
            summary = ", ".join(f"{field}={error:.3e}" for field, error in maxima[group].items())
            print(f"  {group}: {compared}/{attempted} compared | {summary}")

    if failures:
        parser.exit(1, "\n".join(failures) + "\n")


if __name__ == "__main__":
    main()
