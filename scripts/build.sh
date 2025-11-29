#!/bin/zsh
set -e
set -o pipefail

echo "=== 1. GCC Run ==="
cmake --workflow --preset gcc-run

echo "\n=== 2. Intel Run ==="
cmake --workflow --preset intel-run

echo "\n=== 3. NVIDIA Run ==="
cmake --workflow --preset nvidia-run

echo "\nAll workflows completed successfully!"