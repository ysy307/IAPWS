#!/bin/zsh
set -e
set -o pipefail

for dir in build bin lib; do
    if [[ -d $dir ]]; then
        find "$dir" -mindepth 1 -exec rm -rf {} +
    fi
done
