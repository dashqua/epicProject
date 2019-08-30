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

\*---------------------------------------------------------------------------*/

#include "RiemannSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //



void RiemannSolver::computeFlux
(
     surfaceScalarField& rhoFlux,
     surfaceVectorField& rhoUFlux,
     surfaceScalarField& rhoEFlux,
     volVectorField& gradP,
     volTensorField& gradU,
     volVectorField& gradT
)
{
  /* This function compute the 3 fluxes using Riemann Solver.
     It is inspired from numericFlux::computeFlux() from dbnsFlux library
     in foam-extend-4.0. This function does not require foam-extend-4.0 but
     OpenFOAM-6 instead.
  */

  // Step 1a: Settings
  const unallocLabelList& owner          = mesh_.owner();
  const unallocLabelList& neighbour      = mesh_.neighbour();
  const surfaceVectorField& Sf           = mesh_.Sf();
  const surfaceScalarField& magSf        = mesh_.magSf();
  const volVectorField& cellCentre       = mesh_.C();
  const surfaceVectorField& faceCentre   = mesh_.Cf();
  const volScalarField Cv                = thermo_.Cv();
  const volScalarField R                 = thermo_.Cp() - Cv;
  gradP = fvc::grad(p_);
  gradU = fvc::grad(U_);
  gradT = fvc::grad(T_);

  // Step 1b: Limiters stuff

  // Step 2a: Computation of fluxes for internal faces
  forAll(owner, face)
  {
    const label own = owner[face];
    const label nei = owner[face];
    const vector deltaRLeft  = faceCentre[face] - cellCentre[own];
    const vector deltaRRight = faceCentre[face] - cellCentre[nei];

    evaluateFlux
    (

    );
  }

  // Step 2b: Computation of fluxes for boundary faces
  forAll(rhoFlux.boundaryField(), patch)
  {
    const fvPatch& curPatch = p_.boundaryField()[patch].patch();
    // Fluxes
    fvsPatchScalarField& pRhoFlux  = rhoFlux.boundaryField()[patch];
    fvsPatchVectorField& pRhoUFlux = rhoUFlux.boundaryField()[patch];
    fvsPatchScalarField& pRhoEFLux = rhoEFlux.boundaryField()[patch];
    // Patch Fields
    const fvPatchScalarField& pp   = p_.boundaryField()[patch];
    const vectorField& pU          = U_.boundaryField()[patch];
    const scalarField& pT          = T_.boundaryField()[patch];
    const scalarField& pCv         = Cv.boundaryField()[patch];
    const scalarField& pR          = R.boundaryField()[patch];
    //Gradients
    const fvPatchVectorField& pGradP = gradP.boundaryField()[patch];
    const fvPatchTensorField& pGradU = gradU.boundaryField()[patch];
    const fvPatchVectorField& pGradT = gradT.boundaryField()[patch];
    //Limiters stuff
    //Face areas
    const fvsPatchVectorField& pSf = Sf.boundaryField()[patch];
    const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patch];

    forAll(pp, facei)
    {
      evaluateFlux
      (

      );
    }
  }    
}


RiemannSolver::evaluateFlux();
{


}









void RiemannSolver::evaluateFlux(){}









// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::RiemannSolver::operator=(const RiemannSolver& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}
*/
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
