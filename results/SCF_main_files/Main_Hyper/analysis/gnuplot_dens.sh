#!/bin/bash
nam=$(echo $1)
gnuplot -p -e "nam='$nam'" dens.gnu
