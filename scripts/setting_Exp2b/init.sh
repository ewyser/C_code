#!/bin/bash
clear
rm *.txt
rm *.png
rm ep23De_cpu
julia -O3 ./scripts/init.jl
