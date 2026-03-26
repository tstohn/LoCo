Write-Host "=== LoCo Bootstrap ==="

# -----------------------
# 1. Clone vcpkg if missing
# -----------------------
if (!(Test-Path "vcpkg")) {
    Write-Host "Cloning vcpkg..."
    git clone https://github.com/microsoft/vcpkg.git
}

cd vcpkg

# Only unshallow if shallow
if ((git rev-parse --is-shallow-repository) -eq "true") {
    git fetch --unshallow
}
git fetch --all

# Bootstrap vcpkg (disable telemetry)
Write-Host "Bootstrapping vcpkg..."
cmd /c "bootstrap-vcpkg.bat --disableMetrics"

# Install packages using CMD
Write-Host "Installing vcpkg packages..."
cmd /c "vcpkg.exe install"

cd ..

# -----------------------
# 2. Clone nanoflann if missing
# -----------------------
$NanoDir = Join-Path $PWD "dependencies\nanoflann"
if (!(Test-Path $NanoDir)) {
    Write-Host "Cloning nanoflann v1.3.2..."
    git clone https://github.com/jlblancoc/nanoflann.git --branch v1.3.2 $NanoDir
}

# -----------------------
# 3. Build project with CMake
# -----------------------
$NanoInclude = Join-Path $NanoDir "include"

Write-Host "Configuring CMake..."
cmake -S . -B build `
      -DNANOFLANN_INCLUDE_DIR="$NanoInclude" `
      -DCMAKE_TOOLCHAIN_FILE="$PWD/vcpkg/scripts/buildsystems/vcpkg.cmake" `
      -G Ninja `
      -DCMAKE_BUILD_TYPE=Release

Write-Host "Building project..."
cmake --build build --parallel

Write-Host "=== Build Complete ==="