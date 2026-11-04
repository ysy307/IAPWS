"""Shape failures must be deliberate API errors, not bounds faults or crashes."""
import subprocess
import sys

for model in ("97", "95", "ice", "sat"):
    for mode in (("output", "status") if model == "sat" else ("input", "output", "status")):
        result = subprocess.run([sys.argv[1], model, mode], capture_output=True, text=True)
        if result.returncode == 0 or "IAPWS batch:" not in result.stdout + result.stderr:
            raise AssertionError(f"{model}/{mode}: missing shape rejection\n{result.stdout}\n{result.stderr}")
print("All models reject inconsistent shapes before evaluation")
