"""Exercise validation exit status and keep all fixture/log writes in scratch space."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    executable = Path(sys.argv[1]).resolve()
    fixtures = Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="iapws-validation-") as scratch:
        scratch = Path(scratch)
        for case, replacement in (("valid", None), ("mismatch", "999"), ("nan", "NaN"), ("infinity", "Infinity")):
            work = scratch / case
            work.mkdir()
            data = work / "data"
            shutil.copytree(fixtures, data, ignore=shutil.ignore_patterns("*.log"))
            if replacement is not None:
                path = data / "iapws95" / "Helmholtz.dat"
                rows = path.read_text().splitlines()
                fields = rows[1].split()
                fields[2] = replacement
                rows[1] = " ".join(fields)
                path.write_text("\n".join(rows) + "\n")
            result = subprocess.run([str(executable), str(data)], cwd=work, capture_output=True, text=True)
            if (result.returncode == 0) != (replacement is None):
                raise AssertionError(f"{case}: unexpected exit {result.returncode}\n{result.stdout}\n{result.stderr}")
            logs = [work / f"test_iapws{model}.log" for model in ("95", "97", "06")]
            if not all(path.is_file() for path in logs):
                raise AssertionError(f"{case}: logs must be written in the working directory")
            if replacement is not None and "**FAIL**" not in logs[0].read_text():
                raise AssertionError(f"{case}: failure was not recorded")
            if list(data.rglob("*.log")):
                raise AssertionError(f"{case}: validation must not write beside reference data")
    print("Validation harness accepts reference data and rejects mismatch/NaN/infinity")


if __name__ == "__main__":
    main()
