#!/bin/zsh
set -e
set -o pipefail

# Keep the original CMakeBuild.sh untouched.
# This script runs fast release workflows directly.

LOG_DIR="/workspaces/IAPWS/validations"
LOG_FILES=($LOG_DIR/**/*.log(N))
if (( ${#LOG_FILES[@]} > 0 )); then
    rm -f $LOG_FILES
fi

# Use all available cores unless caller overrides it.
: ${CMAKE_BUILD_PARALLEL_LEVEL:=$(nproc)}
export CMAKE_BUILD_PARALLEL_LEVEL

echo "=== 1. Intel (Release/Fast) ==="
cmake --workflow --preset workflow-intel-release

echo "=== 2. GCC (Release/Fast) ==="
cmake --workflow --preset workflow-gcc-release

echo "=== 3. NVIDIA (Release/Fast) ==="
cmake --workflow --preset workflow-nvidia-release

echo "All fast workflows completed successfully!"
