#!/usr/bin/env bash
set -e
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"

echo "=== LoCo Bootstrap ==="

OS="$(uname -s)"
echo "Detected: $OS"

# -------------------------------
# Install base tools
# -------------------------------
if [[ "$OS" == "Linux" ]]; then
    if command -v apt >/dev/null; then
        sudo apt update
        sudo apt install -y git cmake ninja-build build-essential curl zip unzip tar pkg-config
    elif command -v dnf >/dev/null; then
        sudo dnf install -y git cmake ninja-build gcc-c++ make curl zip unzip tar pkgconf-pkg-config
    elif command -v pacman >/dev/null; then
        sudo pacman -S --noconfirm git cmake ninja base-devel curl zip unzip tar pkgconf
    fi
elif [[ "$OS" == "Darwin" ]]; then
    if ! command -v brew >/dev/null; then
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
    brew install git cmake ninja
fi

# -------------------------------
# vcpkg
# -------------------------------
if [ ! -d "vcpkg" ]; then
    echo "Cloning vcpkg..."
    # We clone the full history so baselines can be checked out
    git clone https://github.com/microsoft/vcpkg.git
    cd vcpkg
    ./bootstrap-vcpkg.sh
    cd ..
else
    echo "vcpkg already exists, updating..."
    cd vcpkg
    git fetch --all
    cd ..
fi

# -------------------------------
# Build
# -------------------------------
cmake --preset release
cmake --build --preset release -j

echo "=== DONE ==="
echo "Binary: ./bin/loco"