#!/bin/sh

cp geometry/dom.msh .
gmshToFoam dom.msh
rm dom.msh
paraFoam -builtin -touch
