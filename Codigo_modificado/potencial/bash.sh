#!/bin/bash
#BSUB -oo x.o
#BSUB -eo x.e
date
gfortran -O3 npt-PA10_mod.f -o npt-PA10_mod && ./npt-PA10_mod
cd graficas/
gnuplot elipses_final.gnu
gnuplot elipses_inicio.gnu
gnuplot elipses_intermedia.gnu
gnuplot clusters.gnu
python3 EvsN.py
cd gif/
python3 elipses.py
date