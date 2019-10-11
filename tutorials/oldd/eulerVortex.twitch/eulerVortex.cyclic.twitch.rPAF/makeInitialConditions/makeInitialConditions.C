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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0))
    );
*/

  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 1)
   );

  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(1,1,0))
  );
  
  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 1)
  );

  scalar uInf = 1, vInf = 1, TInf = 1;
  scalar gamma = 1.4, eps = 5;
  scalar pi = constant::mathematical::pi, pi2 = pi*pi;
  scalar x0 = 0, y0 = 0;
  forAll(mesh.cells(), cell)
  {
    scalar x = mesh.C()[cell].x() , y = mesh.C()[cell].y();
    scalar r = Foam::sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
    T[cell] = TInf - (gamma-1)*eps*eps/(8*gamma*pi2) * Foam::exp(1-r*r);
    p[cell] = Foam::pow( T[cell] , gamma/(gamma-1));
    U[cell] = vector(
		     uInf-eps/(2*pi)*Foam::exp((1-r*r)/2)*(y-y0),
		     vInf+eps/(2*pi)*Foam::exp((1-r*r)/2)*(x-x0),
		     0
		     );
  }


  U.write();
  T.write();
  p.write();



    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
