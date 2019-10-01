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
#include "dynamicFvMesh.H"
#include "dynamicArbitraryFvMesh.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createDynamicFvMesh.H"
#include "createDyMControls.H"
  
  Info << "\nCorrection of Initial Variables\n" << endl;
  Info << " [+] Correcting theoretical variables\n";
  Info << " [+] Correcting U, rho and p fields\n";
  scalar M  = 0.5;
  scalar M2 = M*M;
  scalar I  = 5.0;
  scalar I2 = I*I;
  scalar r  = 1.5;
  scalar r2 = r*r;
  scalar theta = Foam::atan(0.5);
  scalar Rhoinf = 1;
  scalar Uinf = 0.8944;
  scalar Vinf = 0.4472;
  scalar Pinf = 3.;
  scalar pii  = constant::mathematical::pi;
  scalar pii2 = pii*pii;
  scalar t = mesh.time().value();
  volVectorField C = mesh.C();

  /*
  const polyBoundaryMesh& bm = mesh.boundaryMesh();
  const label& inletPatchID = bm.findPatchID("inlet");
  const label& outletPatchID = bm.findPatchID("outlet");
  const label& sidesPatchID = bm.findPatchID("emptyPatches");
  const label& topPatchID = bm.findPatchID("top");
  const label& bottomPatchID = bm.findPatchID("bottom");
  */
  


  
  Info << " [+] Instanciation of boundary conditions for rho, p, U, T\n";
  volScalarField rho
  (
   IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("rho", dimensionSet(1,-3,0,0,0,0,0), 1)
   );

  /*
  rho.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[inletPatchID], rho));
  rho.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[outletPatchID], rho));
  rho.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], rho));
  rho.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[topPatchID], rho));
  rho.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[bottomPatchID], rho));  
  */
  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 1)
   );
  /*
  p.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[inletPatchID], p));
  p.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[outletPatchID], p));
  p.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], p));
  p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[topPatchID], p));
  p.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[bottomPatchID], p));  
  */
  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0))
  );
  /*
  U.boundaryFieldRef().set(inletPatchID,fvPatchField<vector>::New("cyclic", mesh.boundary()[inletPatchID], U));
  U.boundaryFieldRef().set(outletPatchID,fvPatchField<vector>::New("cyclic", mesh.boundary()[outletPatchID], U));
  U.boundaryFieldRef().set(sidesPatchID,fvPatchField<vector>::New("cyclic", mesh.boundary()[sidesPatchID], U));
  U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("cyclic", mesh.boundary()[topPatchID], U));
  U.boundaryFieldRef().set(bottomPatchID,fvPatchField<vector>::New("cyclic", mesh.boundary()[bottomPatchID], U));  
  */
  /*
  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 293.15)
  );
  */
  /*
  T.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[inletPatchID], T));
  T.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[outletPatchID], T));
  T.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], T));
  T.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[topPatchID], T));
  T.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("cyclic", mesh.boundary()[bottomPatchID], T));  
  */

  /*
  Info << " [+] loop on patches for specific fixedValue condition\n";

  forAll( bm[inletPatchID] , iface)  // INLET conditions
    {
      scalar y = mesh.Cf().boundaryField()[inletPatchID][iface].y();
      U.boundaryFieldRef()[inletPatchID][iface] = vector(-Foam::exp(-y*y)*Foam::cos(t)*Foam::cos(t+1),0,0);
    }

  forAll( bm[outletPatchID] , iface)  // OUTLET conditions
    {
      scalar y = mesh.Cf().boundaryField()[outletPatchID][iface].y();
      U.boundaryFieldRef()[outletPatchID][iface] = vector(-Foam::exp(-y*y)*Foam::cos(t)*Foam::cos(t+1),0,0);
      p.boundaryFieldRef()[outletPatchID][iface] = 10;
    }

  forAll( bm[sidesPatchID] , iface)  // SIDES conditions
    {
      scalar y = mesh.Cf().boundaryField()[sidesPatchID][iface].y();

    }

  forAll( bm[topPatchID] , iface)  // TOP conditions
    {
      scalar y = mesh.Cf().boundaryField()[topPatchID][iface].y();

    }

  forAll( bm[bottomPatchID] , iface)  // BOTTOM conditions
    {
      scalar y = mesh.Cf().boundaryField()[bottomPatchID][iface].y();

    }
  */


  autoPtr<psiThermo> thermo( psiThermo::New(mesh) );
  volScalarField Cv_ = thermo->Cv();
  volScalarField Cp_ = thermo->Cp();

  forAll(C, cell)
    {
      scalar gamma = Cp_[cell]/Cv_[cell];
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);      

      // Correcting theoretical U, rho and p 
      rho[cell] = Rhoinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , 1/(gamma-1)); 
      U[cell] = vector(
             Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             0);
      p[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , gamma/(gamma-1));
    }


  /*
  U.correctBoundaryConditions();
  p.correctBoundaryConditions();
  T.correctBoundaryConditions();
  */
  volScalarField& T = thermo->T();
  
  U.write();
  p.write();
  T.write();



  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
