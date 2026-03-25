#!/usr/bin/env bash
set -e

echo "=== LoCo Bootstrap ==="

# 1. Fix PATH so tools are found (Homebrew / Linux)
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"

OS="$(uname -s)"
echo "Detected: $OS"

# -------------------------------
# 2. Install base tools
# -------------------------------
if [[ "$OS" == "Linux" ]]; then
    if command -v sudo >/dev/null; then
        sudo apt-get update
        sudo apt-get install -y git cmake ninja-build build-essential curl zip unzip tar pkg-config
    else
        echo "Skipping system install; assuming git, cmake, ninja exist"
    fi
elif [[ "$OS" == "Darwin" ]]; then
    brew install git cmake ninja || echo "Skipping install; tools may already exist"
fi

# -------------------------------
# 3. vcpkg
# -------------------------------
if [ ! -d "vcpkg" ]; then
    git clone https://github.com/microsoft/vcpkg.git
else
    cd vcpkg
    git fetch --unshallow || git fetch --all
    cd ..
fi

cd vcpkg
./bootstrap-vcpkg.sh
cd ..

# -------------------------------
# 4. CMake build
# -------------------------------
cmake --preset release -G Ninja \
      -DCMAKE_TOOLCHAIN_FILE=$(pwd)/vcpkg/scripts/buildsystems/vcpkg.cmake

cmake --build --preset release --parallel

echo "=== DONE ==="
echo "Binary: ./bin/loco"