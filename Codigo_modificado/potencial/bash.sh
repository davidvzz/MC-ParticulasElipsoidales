#!/bin/bash
#BSUB -oo x.o
#BSUB -eo x.e
date
gfortran npt-PA10_mod.f -o npt-PA10_mod && "/home/mejia/Documents/CodigoFortran/Codigo_modificado/potencial/"npt-PA10_mod
cd /home/mejia/Documents/CodigoFortran/Codigo_modificado/potencial/graficas/
gnuplot elipses_final.gnu
gnuplot elipses_inicio.gnu
gnuplot elipses_intermedia.gnu
gnuplot clusters.gnu
python3 EvsN.py
date