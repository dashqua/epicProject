#!/bin/sh

echo "Preparation stage"
echo " [+] Cleaning"
rm *{0123456789}.* log* *.OpenFOAM -rf 
rm constant/poly* -rf

echo " [+] Import basic geometry"
cp geometry/box.msh .
gmshToFoam box.msh > log.gmshToFoam
rm box.msh

echo " [+] Create Paraview touch file"
paraFoam -touch > /dev/null

echo " [+] Initiate Cyclic Patches"
createPatch > log.createPatch
mkdir -p obj
mv *.obj obj/
