Write-Host "=== LoCo Bootstrap ==="

# Clone vcpkg if missing
if (!(Test-Path "vcpkg")) {
    git clone https://github.com/microsoft/vcpkg.git
}

cd vcpkg

# Only unshallow if shallow
if ((git rev-parse --is-shallow-repository) -eq "true") {
    git fetch --unshallow
}
git fetch --all

# Bootstrap vcpkg (disable telemetry)
cmd /c "bootstrap-vcpkg.bat --disableMetrics"

# Install packages using CMD, not PowerShell parsing
cmd /c "vcpkg.exe install"

cd ..