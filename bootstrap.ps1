Write-Host "=== LoCo Bootstrap ==="

# Clone vcpkg if missing
if (!(Test-Path "vcpkg")) {
    git clone https://github.com/microsoft/vcpkg.git
}

cd vcpkg

# Make sure we have the full history
git fetch --unshallow
git fetch --all

# Bootstrap vcpkg (disable telemetry)
.\bootstrap-vcpkg.bat --disableMetrics

#install dependencies from manifest
.\vcpkg.exe install
cd ..

# Configure & build with Visual Studio generator
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=.\vcpkg\scripts\buildsystems\vcpkg.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release

Write-Host "=== DONE ==="
Write-Host "Binary: .\build\Release\loco.exe"