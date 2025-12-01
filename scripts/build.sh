#!/bin/zsh
set -e
set -o pipefail

LOG_DIR="/workspaces/IAPWS/validations"
LOG_FILES=($LOG_DIR/**/*.log(N))
if (( ${#LOG_FILES[@]} > 0 )); then
    rm -f $LOG_FILES
fi

echo "=== 1. Intel (Debug) ==="
cmake --workflow --preset workflow-intel-debug

echo "=== 2. GCC (Debug) ==="
cmake --workflow --preset workflow-gcc-debug

echo "=== 3. NVIDIA (Debug) ==="
cmake --workflow --preset workflow-nvidia-debug

echo "All workflows completed successfully!"
