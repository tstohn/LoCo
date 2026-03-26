#!/usr/bin/env bash
set -e

echo "=== LoCo Bootstrap ==="

OS="$(uname -s)"
echo "Detected OS: $OS"

# -----------------------
# 1. Install base tools
# -----------------------
if [[ "$OS" == "Linux" ]]; then
    if command -v sudo >/dev/null; then
        sudo apt-get update
        sudo apt-get install -y build-essential cmake ninja-build gcc gfortran git curl unzip zip pkg-config
    else
        echo "Skipping system install; ensure git, cmake, ninja, compiler exist"
    fi
elif [[ "$OS" == "Darwin" ]]; then
    brew install gcc gfortran git cmake ninja || echo "Tools may already exist"
    export PATH="/opt/homebrew/bin:$PATH"   # Apple Silicon
    export PATH="/usr/local/bin:$PATH"      # Intel mac
    export PATH="/opt/homebrew/opt/gcc/bin:$PATH"
fi

# -----------------------
# 2. vcpkg clone / bootstrap
# -----------------------
if [ ! -d "vcpkg" ]; then
    git clone https://github.com/microsoft/vcpkg.git
else
    cd vcpkg
    git fetch --unshallow || git fetch --all
    cd ..
fi

cd vcpkg
git fetch --all
git checkout master
./bootstrap-vcpkg.sh
cd ..

# -----------------------
# 3. Build project
# -----------------------
cmake -S . -B build \
      -DCMAKE_TOOLCHAIN_FILE=$(pwd)/vcpkg/scripts/buildsystems/vcpkg.cmake \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build --parallel

echo "=== Build Complete ==="