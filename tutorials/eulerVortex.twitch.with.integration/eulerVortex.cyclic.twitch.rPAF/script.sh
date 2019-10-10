#!/usr/bin/sh

echo "\nPreparation stage"
echo " [+] clean the case"
rm -rf  ./0*/ ./1/

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


