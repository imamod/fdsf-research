#! /bin/bash

BUILD_TYPE=$1

mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE
cmake --build .
