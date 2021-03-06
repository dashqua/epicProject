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

  scalar t = mesh.time().value();
  volVectorField C = mesh.C();



  
  
  const polyBoundaryMesh& bm = mesh.boundaryMesh();
  const label& inletPatchID = bm.findPatchID("inlet");
  const label& outletPatchID = bm.findPatchID("outlet");
  const label& sidesPatchID = bm.findPatchID("sides");
  const label& topPatchID = bm.findPatchID("top");
  const label& bottomPatchID = bm.findPatchID("bottom");

  

  Info << " [+] Instanciation of boundary conditions for p, U, T\n";


    
  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 1)
   );
  
  p.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[inletPatchID], p));
  p.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[outletPatchID], p));
  p.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], p));
  p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[topPatchID], p));
  p.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[bottomPatchID], p));  



  
  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0))
  );
  U.boundaryFieldRef().set(inletPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[inletPatchID], U));
  U.boundaryFieldRef().set(outletPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[outletPatchID], U));
  U.boundaryFieldRef().set(sidesPatchID,fvPatchField<vector>::New("empty", mesh.boundary()[sidesPatchID], U));
  U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("zeroGradient", mesh.boundary()[topPatchID], U));
  U.boundaryFieldRef().set(bottomPatchID,fvPatchField<vector>::New("zeroGradient", mesh.boundary()[bottomPatchID], U));  


  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 293.15)
  );
  T.boundaryFieldRef().set(inletPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[inletPatchID], T));
  T.boundaryFieldRef().set(outletPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[outletPatchID], T));
  T.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], T));
  T.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[topPatchID], T));
  T.boundaryFieldRef().set(bottomPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[bottomPatchID], T));  

  
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



  Info << " [+] Set up of internal fields\n";


  // Internal Field
  forAll(mesh.C(), cell)   // loop over cell
    {
      U[cell] = vector(0,0,0);
      p[cell] = 0;
      T[cell] = 293.15;
    }


  Info << " [+] Write-up of variables\n";


  U.write();
  p.write();
  T.write();


  
  Info << " [+] Correction of variables using thermo model\n";
  

  autoPtr<fluidThermo> pThermo( fluidThermo::New(mesh) );
  fluidThermo& thermo = pThermo();
  //thermo.validate(args.executable(), "h", "e");
  //volScalarField& pT = thermo.p();
  //pT.write();
  //T.write();
  
  

  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
