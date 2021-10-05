cmake_minimum_required(VERSION 3.0)

# set compiler and project name
set(CMAKE_C_COMPILER "gcc")
project(ep23De_cpu VERSION 1.0 LANGUAGES C)

add_compile_options(-w -O3 -DDAT=double -DSIM=1 -DFPS=5 -DD=0.1)
add_subdirectory(functions)
add_executable(ep23De_cpu main.c)
target_link_libraries(ep23De_cpu PUBLIC functions)



















#https://cmake.org/cmake/help/latest/guide/tutorial/Adding%20a%20Library.html




