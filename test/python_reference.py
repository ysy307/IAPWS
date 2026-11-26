"""Compare the Fortran kernels with a live iapws 1.5.5 calculation."""
import math
import subprocess
import sys

try:
    from iapws import IAPWS95, IAPWS97
    from iapws import _iapws, iapws95
    from scipy.optimize import fsolve, newton
    from iapws._iapws import _Ice, _Supercooled
    from iapws.iapws97 import _Region2_meta
except ImportError:
    sys.exit(77)


def reference(model, temperature, pressure):
    p = pressure * 1e-6
    if model == "sat":
        ref = IAPWS95(T=temperature, x=0.5)
        return {0: ref.P * 1e6, 1: ref.Liquid.rho, 2: ref.Vapor.rho,
                3: ref.Liquid.h * 1e3, 4: ref.Vapor.h * 1e3,
                5: ref.Liquid.s * 1e3, 6: ref.Vapor.s * 1e3}
    if model == "g12":
        ref = _Supercooled(temperature, p)
        return {0: ref["rho"], 3: ref["cp"] * 1e3, 4: ref["w"],
                7: ref["alfap"], 8: ref["xkappa"] * 1e-6}
    if model == "g9":
        # The residual terms are negligible at this density.
        tau = 647.096 / temperature
        ideal = IAPWS95()._phi0(tau, pressure / (461.51805 * temperature * 322))
        return {3: 461.51805 * (1 - tau**2 * ideal["fiott"])}
    if model == "ice":
        ref = _Ice(temperature, p)
        return {0: ref["rho"], 1: ref["h"] * 1e3, 2: ref["s"] * 1e3, 3: ref["cp"] * 1e3}
    if model == "meta":
        ref = _Region2_meta(temperature, p)
        return {0: 1 / ref["v"], 1: ref["h"] * 1e3, 2: ref["s"] * 1e3,
                3: ref["cp"] * 1e3, 4: ref["w"]}
    ref = (IAPWS95 if model == "95" else IAPWS97)(T=temperature, P=p)
    if model == "95":
        # Match the native (T,rho) EOS evaluation after finding the reference
        # density. The Python TP object also uses prescribed pressure in some
        # derivative expressions and would otherwise amplify its solve residual.
        ref = IAPWS95(T=temperature, rho=ref.rho)
    return {0: ref.rho, 1: ref.h * 1e3, 2: ref.s * 1e3, 3: ref.cp * 1e3,
            4: ref.w, 5: ref.mu, 6: ref.k}


# Refine implicit reference solves to the numerical floor instead of accepting
# SciPy's default, which is suited to ordinary engineering precision.
def ref_newton(*args, **kwargs):
    kwargs.setdefault("tol", 4e-16)
    return newton(*args, **kwargs)


def ref_fsolve(*args, **kwargs):
    kwargs.setdefault("xtol", 5e-14)
    return fsolve(*args, **kwargs)


_iapws.newton = ref_newton
iapws95.fsolve = ref_fsolve

cases = [
    ("95", 300.0, 1e5), ("95", 400.0, 1e5), ("95", 650.0, 25e6),
    ("97", 300.0, 3e6), ("97", 700.0, 3e6), ("97", 650.0, 25e6),
    ("97", 1500.0, 30e6),
    ("ice", 260.0, 1e5), ("meta", 450.0, 1e6),
    ("g12", 273.15, 101325.0), ("g12", 235.15, 101325.0),
    ("g12", 250.0, 200e6), ("g12", 200.0, 400e6), ("g12", 250.0, 400e6),
    ("g9", 50.0, 1e-42), ("g9", 100.0, 1e-16),
    ("sat", 275.0, 0.0), ("sat", 450.0, 0.0), ("sat", 625.0, 0.0),
]
data = "".join(f"{model} {temperature} {pressure}\n" for model, temperature, pressure in cases)
process = subprocess.run([sys.argv[1]], input=data, text=True, capture_output=True)
if process.returncode:
    raise AssertionError(process.stdout + process.stderr)
lines = process.stdout.splitlines()
if len(lines) != len(cases):
    raise AssertionError(f"expected {len(cases)} states, got {len(lines)}")

max_error = {}
for (model, temperature, pressure), line in zip(cases, lines):
    values = list(map(float, line.split()))
    for i, expected in reference(model, temperature, pressure).items():
        if model == "97" and temperature > 1173.15 and i in (5, 6):
            if not math.isnan(values[i]):
                raise AssertionError("transport must reject states outside the official range")
            continue
        relative = abs(values[i] / expected - 1)
        max_error[model] = max(max_error.get(model, 0), relative)
        # The reference EOS sum, implicit density and coexistence solves, and
        # near-zero enthalpy/entropy amplify roundoff differently by state.
        tolerance = {"95": 4e-12, "97": 3e-12, "ice": 3e-13, "meta": 3e-13,
                     "g12": 5e-13, "g9": 5e-14, "sat": 8e-12}[model]
        if not math.isfinite(relative) or relative > tolerance:
            raise AssertionError(
                f"{model} T={temperature} p={pressure} field={i}: "
                f"Fortran={values[i]:.17g}, Python={expected:.17g}, relative={relative:.3g}"
            )
for model, relative in max_error.items():
    print(f"{model}: maximum relative error {relative:.3g}")
