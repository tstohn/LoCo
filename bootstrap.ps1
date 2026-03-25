Write-Host "=== LoCo Bootstrap ==="

if (!(Test-Path "vcpkg")) {
    git clone https://github.com/microsoft/vcpkg.git
    cd vcpkg
    .\bootstrap-vcpkg.bat
    cd ..
}

cmake --preset release
cmake --build --preset release

Write-Host "=== DONE ==="
Write-Host "Binary: .\bin\loco.exe"