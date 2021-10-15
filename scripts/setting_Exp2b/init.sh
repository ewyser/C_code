#!/bin/bash
clear
path="/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code/scripts/setting_Exp2b"
cd ${path}
rm *.txt
rm *.png
rm ep23De_cpu
julia -O3 ./scripts/init.jl
