#!/bin/zsh

ROOT_DIR=$(realpath .)

THIRD_PARTY_DIR="$ROOT_DIR/third_party"
FPM_DIR="$THIRD_PARTY_DIR/fpm"
INSTALL_PREFIX="$THIRD_PARTY_DIR/.local"

if [ ! -d "$THIRD_PARTY_DIR" ]; then
    mkdir -p "$THIRD_PARTY_DIR"
fi
cd "$THIRD_PARTY_DIR"

if [ -d "$FPM_DIR" ]; then
    rm -rf "$FPM_DIR"
fi

git clone https://github.com/fortran-lang/fpm
cd "$FPM_DIR"

./install.sh --prefix="$INSTALL_PREFIX"

export PATH="$INSTALL_PREFIX/bin:$PATH"