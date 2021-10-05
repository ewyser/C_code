#!/bin/bash
clear
#------------------------------------------------------------------
# Refs.
# https://cmake.org/cmake/help/latest/manual/cmake.1.html
#------------------------------------------------------------------
path="/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code"
cd ${path}
echo "============================================================="
echo "CMake: "${path}
echo "-------------------------------------------------------------"
#------------------------------------------------------------------
# remove what should be removed
rm -r build
mkdir ./build

# generate a project buildsystem
cmake -S ./src -B ./build -DCMAKE_BUILD_TYPE=Release
# buil a project
cmake --build ./build

# copy exec code to script
path="./scripts/setting_Exp2b/" 
cp -R ./build/ep23De_cpu ${path}ep23De_cpu
# execute code
cd ${path}
./ep23De_cpu
# display results calling python script
python3 display.py