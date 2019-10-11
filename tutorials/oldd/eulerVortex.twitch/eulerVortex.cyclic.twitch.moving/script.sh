#!/usr/bin/sh

echo "\nPreparation stage"
echo " [+] clean the case"
rm -rf  ./0*/ ./1*/ ./2*/ ./3*/ ./4*/ ./5*/

echo " [+] import geometry"
rm -rf log *.foam 0* constant/polyMesh
mkdir -p log
cp geometry/box.msh .
gmshToFoam  box.msh > log/log.gmshToFoam
rm box.msh
paraFoam -touch -builtin > /dev/null

echo " [+] set empty patches and cyclic boundaries"
createPatch > log/log.createPatch
rm constant/polyMesh -rf
mv 0.001/polyMesh constant/
rm 0.001 -rf

echo " [+] initialize BCs and internal fields"
./makeInitialConditions/makeInitialConditions > log/log.makeInitialConditions


echo "\nSimulation stage"
echo " [+] Run rhoCentralArbitraryMFoam"
rhoCentralArbitraryMFoam > log/log.rhoCentralFoam

#echo " [+] Post-processing"
#gnuplot -p -e "set term png; \
#               set output 'postProcessing/p.png'; \
#               set style line 1 lw 4 lc rgb '#990042' ps 2 pt 6 pi 5; \
#               p 'postProcessing/probes/0/p' u 1:2"

