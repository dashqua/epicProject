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
1;5202;0c
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

  
  Info << " [+] Instanciation of boundary conditions for rho, p, U, T\n";
  volScalarField rho
  (
   IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("rho", dimensionSet(1,-3,0,0,0,0,0), 1)
   );

  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0)
   );

  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0))
  );
  
  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 300)
  );
  
  
  const polyBoundaryMesh& bm = mesh.boundaryMesh();
  const label& sidesPatchID = bm.findPatchID("emptyPatches_empt");
  U.boundaryFieldRef().set(sidesPatchID,fvPatchField<vector>::New("empty", mesh.boundary()[sidesPatchID], U));
  p.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], p));
  T.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], T));
  rho.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], rho));

  const label& topPatchID = bm.findPatchID("top_cyc");
  U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("noSlip", mesh.boundary()[topPatchID], U));
  p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[topPatchID], p));
  T.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[topPatchID], T));
  rho.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[topPatchID], rho));

  const label& bottomPatchID = bm.findPatchID("bottom_cyc");
  U.boundaryFieldRef().set(bottomPatchID,fvPatchField<vector>::New("noSlip", mesh.boundary()[bottomPatchID], U));
  p.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[bottomPatchID], p));
  T.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[bottomPatchID], T));
  rho.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[bottomPatchID], rho));

  const label& inletPatchID = bm.findPatchID("inlet_cyc");
  U.boundaryFieldRef().set(inletPatchID,fvPatchField<vector>::New("noSlip", mesh.boundary()[inletPatchID], U));
  p.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[inletPatchID], p));
  T.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[inletPatchID], T));
  rho.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[inletPatchID], rho));

  const label& outletPatchID = bm.findPatchID("outlet_cyc");
  U.boundaryFieldRef().set(outletPatchID,fvPatchField<vector>::New("noSlip", mesh.boundary()[outletPatchID], U));
  p.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[outletPatchID], p));
  T.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[outletPatchID], T));
  rho.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[outletPatchID], rho));


  T.write();

  autoPtr<psiThermo> pThermo
    (
     psiThermo::New(mesh)
     );
  psiThermo& thermo = pThermo();
  volScalarField Cv_ = thermo.Cv();
  volScalarField Cp_ = thermo.Cp();
  /*
  autoPtr<psiThermo> thermo( psiThermo::New(mesh) );
  volScalarField Cv_ = thermo->Cv();
  volScalarField Cp_ = thermo->Cp();
  */
  forAll(C, cell)
    {
      scalar gamma = Cp_[cell]/Cv_[cell];
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);      

      // Correcting theoretical U, rho and p
      T[cell] = Rhoinf *( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) ); 
      rho[cell] = Foam::pow( T[cell] , 1/(gamma-1) ); 
      U[cell] = vector(
             Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             0);
      p[cell] = Foam::pow( rho[cell] , gamma );
    }

  //volScalarField&
  //T = thermo->T();
  
  U.write();
  p.write();
  rho.write();
  T.write();



  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
