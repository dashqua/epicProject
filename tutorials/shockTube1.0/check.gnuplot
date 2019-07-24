set terminal pngcairo
set output "check.png"
set title "p[20] evolution in time"
p 'evolp_dbnsFoam' u 1:2, 'evolp_epicSolver' u 1:2
