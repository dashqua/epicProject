#!/bin/sh

echo "\nPreparation stage"
echo " [+] Cleaning"
rm ./0*/ log* *.OpenFOAM obj/ log/ -rfd 
rm constant/poly* -rf
mkdir -p log

echo " [+] Import basic geometry"
cp geometry/box.msh .
gmshToFoam box.msh > log/log.gmshToFoam
rm box.msh

echo " [+] Initiate empty and cyclic patches"
createPatch > log/log.createPatch
mkdir -p obj
mv *.obj obj/
rm constant/poly* -rf
mv ./0*?/polyMesh constant/
rm ./0*/ -rf
rm *~ -rf

echo " [+] Run makeInitialConditions"
./makeInitialConditions/makeInitialConditions > log/log.makeInitialConditions 

echo "\nSimulation stage"
#echo " [+] Generate Analytic Sol"
#./makeAnalyticSolution/makeAnalyticSolution > log/log.makeAnalyticSolution

echo " [+] Run Solver"
rhoPimpleArbitraryFoam > log.solver
