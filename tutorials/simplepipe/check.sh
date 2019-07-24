#!/bin/bash
epicProject
./Allclean
dbnsFoam
python3 evol_p.py dbnsFoam
./Allclean
epicSolver
python3 evol_p.py epicSolver
gnuplot -p -e "set terminal png; set output 'check.png'; set title 'p[20] evolution in time'; p 'evolp_dbnsFoam' u 1:2, 'evolp_epicSolver' u 1:2"
rm -rf evolp*
