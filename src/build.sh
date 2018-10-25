#! /bin/bash

BUILD_TYPE=Release

# Сборка Utils
cd Utils
bash ./build.sh $BUILD_TYPE

# Сборка fdsf
cd ../fdsf
bash ./build.sh $BUILD_TYPE
echo `pwd`
# Сборка fd-half
cd ../fd-half
python3 build.py
echo `pwd`
# Сборка fd-Jmhalf
cd ../fd-Jmhalf
python3 build.py
