:@echo off

set BUILD_TYPE="%~1"

if [%BUILD_TYPE%] == [] (
    echo Usage:
    echo    build.bat [Debug/Release/DebugMP/ReleaseMP]
    exit 1
)

mkdir build
cd build
cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_BUILD_TYPE=%BUILD_TYPE%
cmake --build .

:@echo on
