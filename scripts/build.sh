#!/bin/zsh
set -e
set -o pipefail

LOG_DIR="/workspaces/IAPWS/validations"
LOG_FILES=($LOG_DIR/**/*.log(N))
if (( ${#LOG_FILES[@]} > 0 )); then
    rm -f $LOG_FILES
fi

echo "=== 1. Intel Run ==="
cmake --workflow --preset intel-run

echo "=== 2. GCC Run ==="
cmake --workflow --preset gcc-run

echo "=== 3. NVIDIA Run ==="
cmake --workflow --preset nvidia-run

echo "All workflows completed successfully!"
