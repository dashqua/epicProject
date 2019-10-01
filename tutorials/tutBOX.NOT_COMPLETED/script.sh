#!/bin/sh

echo "\nPreparation stage"
echo " [+] Cleaning"
rm ./0*/ ./*e-*  log* *.OpenFOAM obj/ log/ constant/poly* -rfd 
mkdir -p log

echo " [+] Import basic geometry"
cp geometry/box.msh .
gmshToFoam box.msh > log/log.gmshToFoam
rm box.msh

echo " [+] Create patches"
#createPatch > log/log.createPatch
#mkdir -p obj
#mv *.obj obj/ -f
#rm constant/poly* -rf
#mv ./0*?/polyMesh constant/
#rm ./0*/ -rf
#rm *~ -rf

echo " [+] Run checkMesh"
checkMesh > log/log.checkMesh

#echo " [+] Run makeInitialConditions"
#./makeInitialConditions/makeInitialConditions > log/log.makeInitialConditions 

echo "\nSimulation stage"
#rhoCentralFoam
#echo " [+] Generate Analytic Sol"
#./makeAnalyticSolution/makeAnalyticSolution > log/log.makeAnalyticSolution

#echo " [+] Run Solver"
#rhoPimpleFoam > log.solver #rhoPimpleArbitraryFoam > log.solver
