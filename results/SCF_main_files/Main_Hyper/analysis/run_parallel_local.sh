#!/bin/bash

threads=$(grep 'threads' input_general | awk '{print $1}')
export OMP_NUM_THREADS=$threads
ulimit -s unlimited
rm -r OUTPUT*

./analysis.exe
