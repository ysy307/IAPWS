#!/bin/bash

# ヘルプ表示関数
usage() {
    echo "Usage: $0 [compiler|clean] [build_type]"
    echo ""
    echo "Commands:"
    echo "  clean    ビルド成果物を削除 (fpm clean + bin/lib削除)"
    echo ""
    echo "Compilers:"
    echo "  gcc      (gfortran)"
    echo "  intel    (ifx)"
    echo "  nvidia   (nvfortran)"
    echo ""
    echo "Build Types:"
    echo "  debug"
    echo "  release"
    echo ""
    echo "Example:"
    echo "  $0 gcc debug"
    echo "  $0 clean"
    exit 1
}

# クリーンアップ処理
if [ "$1" == "clean" ]; then
    echo "Cleaning build artifacts..."
    set -x
    # fpmのビルドディレクトリを削除 (--all でキャッシュ等も含む)
    fpm clean --all
    # CMakeなどで生成された可能性のある bin, lib ディレクトリも削除
    rm -rf bin lib
    set +x
    echo "Done."
    exit 0
fi

# 引数が足りない場合はヘルプを表示
if [ $# -lt 2 ]; then
    usage
fi

COMPILER_INPUT=$1
BUILD_TYPE_INPUT=$2

# 入力を小文字に変換
COMPILER=$(echo "$COMPILER_INPUT" | tr '[:upper:]' '[:lower:]')
BUILD_TYPE=$(echo "$BUILD_TYPE_INPUT" | tr '[:upper:]' '[:lower:]')

# 変数の初期化
FC=""
FLAGS=""

# コンパイラとビルドタイプに応じたフラグ設定
# CMakePresets.json の cacheVariables に対応しています
case "$COMPILER" in
    gcc)
        FC="gfortran"
        # base: -Wall -Wextra -fopenmp
        BASE_FLAGS="-Wall -Wextra -fopenmp"
        
        if [ "$BUILD_TYPE" == "debug" ]; then
            # gcc-debug
            FLAGS="$BASE_FLAGS -g -O0 -fcheck=all -fbacktrace"
        elif [ "$BUILD_TYPE" == "release" ]; then
            # gcc-release
            FLAGS="$BASE_FLAGS -O3 -march=native -ffast-math"
        else
            echo "Error: Unknown build type '$BUILD_TYPE'"
            usage
        fi
        ;;

    intel)
        FC="ifx"
        # base: -fpp -traceback -warn all -qopenmp
        BASE_FLAGS="-fpp -traceback -warn all -qopenmp"

        if [ "$BUILD_TYPE" == "debug" ]; then
            # intel-debug
            FLAGS="$BASE_FLAGS -g -O0 -check all -fpe0"
        elif [ "$BUILD_TYPE" == "release" ]; then
            # intel-release
            FLAGS="$BASE_FLAGS -O3 -xHost"
        else
            echo "Error: Unknown build type '$BUILD_TYPE'"
            usage
        fi
        ;;

    nvidia)
        FC="nvfortran"
        # base: -Minfo=all -mp
        BASE_FLAGS="-Minfo=all -mp"

        if [ "$BUILD_TYPE" == "debug" ]; then
            # nvidia-debug
            FLAGS="$BASE_FLAGS -g -O0 -Mbounds -Mchkptr -traceback -Ktrap=fp"
        elif [ "$BUILD_TYPE" == "release" ]; then
            # nvidia-release
            FLAGS="$BASE_FLAGS -O3 -fast"
        else
            echo "Error: Unknown build type '$BUILD_TYPE'"
            usage
        fi
        ;;

    *)
        echo "Error: Unknown compiler '$COMPILER'"
        usage
        ;;
esac

# 設定内容の表示
echo "========================================================="
echo " FPM Build Configuration"
echo "========================================================="
echo " Compiler   : $FC"
echo " Build Type : $BUILD_TYPE"
echo " Flags      : $FLAGS"
echo "========================================================="

# 実行コマンドの表示と実行 (set -x)
set -x

# ビルド実行
# CMakeのワークフローと同様に Build -> Test の順で実行します
fpm build --compiler "$FC" --flag "$FLAGS"
BUILD_STATUS=$?

if [ $BUILD_STATUS -eq 0 ]; then
    set +x
    echo ""
    echo "Build successful. Running tests..."
    set -x
    fpm test --compiler "$FC" --flag "$FLAGS"
else
    set +x
    echo "Build failed."
    exit $BUILD_STATUS
fi