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
cmake -S ./src -B ./build
# buil a project
cmake --build ./build
#./build/prog
#------------------------------------------------------------------
