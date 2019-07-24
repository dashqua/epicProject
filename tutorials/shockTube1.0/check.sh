#!/bin/bash
epicProject
./Allclean
dbnsFoam
python3 evol_p.py dbnsFoam
./Allclean
epicSolver
python3 evol_p.py epicSolver
gnuplot check.gnuplot
