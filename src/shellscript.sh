#!/bin/bash
clear
#------------------------------------------------------------------
# Refs.
# https://cmake.org/cmake/help/latest/manual/cmake.1.html
#------------------------------------------------------------------
gcc -DDAT=double -DSIM=1 -DFPS=5 -DD=0.1 -O3 -o cpu main.c ./functions/init.c ./functions/saveData.c ./functions/CFL.c ./functions/getG.c ./functions/topol.c ./functions/NdN.c ./functions/basis.c ./functions/accum.c ./functions/solve.c ./functions/FLIP.c ./functions/DM_BC.c ./functions/strains.c ./functions/elast.c ./functions/DPPlast.c ./functions/volLock.c