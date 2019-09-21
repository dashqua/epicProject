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
#include "psiThermo.H"
#include "bound.H"

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

  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 1)
   );

  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(0,0,0,1,0,0,0), 293.15)
   );


  
  Info << "\nCorrection of Initial Variables\n" << endl;
  scalar t = mesh.time().value();
  volVectorField C = mesh.C();


  const unallocLabelList& owner = mesh.owner();  
  forAll( mesh.boundary() , ipatch)  // loop over all boundary patches
    {     
      if (mesh.boundary()[ipatch].name() == "inlet")  // get the whole "inlet" boundary patch
	{
	  forAll(mesh.boundary()[ipatch], iface) // loop over faces of that patch
	    {
	      const label own = owner[iface]; //change U value for owner cell of that face
	      scalar x = mesh.C()[own].x(), y = mesh.C()[own].y();
	      U[own] = vector(1e4*Foam::exp(-(x+y)*(x+y))*Foam::cos(t)*Foam::cos(t),0,0);
	      Info << x << " " << y << " " << U[own].x() << endl;
	    }
	}
    }

  U.correctBoundaryConditions();
  U.write();



  forAll(mesh.C(), cell)   // loop over cell
    {
      p[cell] = 1;
      T[cell] = 293.15;
    }
    p.correctBoundaryConditions();
    p.write();
    T.correctBoundaryConditions();
    T.write();



  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
