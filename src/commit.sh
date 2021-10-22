#!/bin/bash
clear
#------------------------------------------------------------------
# Refs.
# https://cmake.org/cmake/help/latest/manual/cmake.1.html
#------------------------------------------------------------------
path="/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code"
DATE=$(date)
cd ${path}
echo "Directory: "${path} 
echo ${DATE}
echo "Enter message for commit:" 
read commit_message
git add -A
git commit -m "${commit_message}: on ${DATE}"
git push -u origin bspline