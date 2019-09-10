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

  autoPtr<psiThermo> thermo( psiThermo::New(mesh) );
  
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
  volScalarField rho
  (
   IOobject("rho", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("rho", dimensionSet(1,-3,0,0,0,0,0), 0)
  );
  volScalarField Cv_ = thermo->Cv();
  volScalarField Cp_ = thermo->Cp();

  
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
  Info << " [+] Correcting theoretical and conserved variables\n";
  forAll(C, cell)
    {
      scalar gamma = Cp_[cell]/Cv_[cell];
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      //scalar x3 = C[cell].z();
      //vector pos = vector(x1, x2, x3);
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);      
      // Correcting theoretical U, rho and p 
      rho[cell] = Rhoinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , 1/(gamma-1)); Info << "1";
      U[cell] = vector(
             Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             0);
      p[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , gamma/(gamma-1));

      // Correcting Primitive Fields
      //rho_[cell] = rho_theo_[cell];
      //U_[cell]   = U_theo_[cell];
      //p_[cell]   = p_theo_[cell];

      // Correcting Conserved Fields
      //rhoU_[cell] = rho_[cell] * U_[cell];
      //rhoE_[cell] = rho_[cell] * E_[cell]; // (h_[cell] + 0.5*magSqr(U_));
      // Fluxes don't need to be initialized because they are computed afterwards   
    }

  p.correctBoundaryConditions();
  U.correctBoundaryConditions();
  rho.correctBoundaryConditions();
  T.correctBoundaryConditions();
  
  p.write();
  U.write();
  T.write();
  rho.write();






  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
