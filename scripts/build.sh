#!/bin/zsh
set -e
set -o pipefail

LOG_DIR="/workspaces/IAPWS/validations"
# 配列に展開、存在しないファイルは無視
LOG_FILES=($LOG_DIR/**/*.log(N))

if (( ${#LOG_FILES[@]} > 0 )); then
    rm -f $LOG_FILES
fi

echo "=== 1. GCC Run ==="
cmake --workflow --preset gcc-run

echo "\n=== 2. Intel Run ==="
cmake --workflow --preset intel-run

echo "\n=== 3. NVIDIA Run ==="
cmake --workflow --preset nvidia-run

echo "All workflows completed successfully!"
