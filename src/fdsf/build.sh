#! /bin/bash

BUILD_TYPE=$1

if [ -z "$BUILD_TYPE"]; then
    BUILD_TYPE=Release
fi

mkdir -p build/$BUILD_TYPE
cd build
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE
cmake --build .
