/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    makeInitialConditions

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
  
  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0))
  );
  volScalarField T
  (
   IOobject("T", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 273.15)
  );
  volScalarField p
  (
   IOobject("p", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0)
  );
  
  forAll(mesh.C(), cell)
  {
    if (mesh.C()[cell].x() < 0.5)
      p[cell] = 303975;
    if (mesh.C()[cell].x() == 0.5)
      p[cell] = (303975+101325)/2;
    if (mesh.C()[cell].x() > 0.5)
      p[cell] = 101325;

    U[cell] = vector(0,0,0);
    T[cell] = 273.15;
  }

  p.correctBoundaryConditions();
  U.correctBoundaryConditions();
  T.correctBoundaryConditions();
  
  p.write();
  U.write();
  T.write();







  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
