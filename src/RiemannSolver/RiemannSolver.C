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
     It is inspired from numericFlux::computeFlux() of dbnsFlux library
     in foam-extend-4.0. This function should not require foam-extend-4.0 but
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
  //gradP.correctBoundaryConditions();
  //gradU.correctBoundaryConditions();
  //gradT.correctBoundaryConditions();

  // Step 1b: Limiters stuff
  MDLimiter<scalar, BarthJespersenLimiter> scalarPLimiter(p_, gradP);
  MDLimiter<vector, BarthJespersenLimiter> vectorULimiter(U_, gradU);
  MDLimiter<scalar, BarthJespersenLimiter> scalarTLimiter(T_, gradT);
  const volScalarField& pLimiter = scalarPLimiter.phiLimiter();
  const volVectorField& ULimiter = vectorULimiter.phiLimiter();
  const volScalarField& TLimiter = scalarTLimiter.phiLimiter();
  
  
  // Step 2a: Computation of fluxes for internal faces
  forAll(owner, face)
  {
    const label own = owner[face];
    const label nei = owner[face];
    const vector deltaRLeft  = faceCentre[face] - cellCentre[own];
    const vector deltaRRight = faceCentre[face] - cellCentre[nei];
    
    evaluateFluxInternal
    (
     rhoFlux[face],
     rhoUFlux[face],
     rhoEFlux[face],
     p_[own] + pLimiter[own] * (deltaRLeft & gradP[own]),
     p_[nei] + pLimiter[nei] * (deltaRRight & gradP[nei]),
     U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
     U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
     T_[own] + TLimiter[own] * (deltaRLeft & gradT[own]),
     T_[nei] + TLimiter[nei] * (deltaRRight & gradT[nei]),
     R[own],
     R[nei],
     Cv[own],
     Cv[nei],
     Sf[face],
     magSf[face]
    );
    
  }

  // Step 2b: Computation of fluxes for boundary faces
  forAll(rhoFlux.boundaryField(), patch)
  {
    const fvPatch& curPatch = p_.boundaryField()[patch].patch();
    // Fluxes
    fvsPatchScalarField& pRhoFlux  = rhoFlux.boundaryField()[patch];
    fvsPatchVectorField& pRhoUFlux = rhoUFlux.boundaryField()[patch];
    fvsPatchScalarField& pRhoEFlux = rhoEFlux.boundaryField()[patch];
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
    const fvPatchScalarField& pPatchLim = pLimiter.boundaryField()[patch];
    const fvPatchVectorField& UPatchLim = ULimiter.boundaryField()[patch];
    const fvPatchScalarField& TPatchLim = TLimiter.boundaryField()[patch];    
    //Face areas
    const fvsPatchVectorField& pSf = Sf.boundaryField()[patch];
    const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patch];
    //
    forAll(pp, face)
    {
      evaluateFluxBoundary
      (
       pRhoFlux[face],
       pRhoUFlux[face],
       pRhoEFlux[face],
       pp[face],
       pp[face],
       pU[face],
       pU[face],
       pT[face],
       pT[face],
       pR[face],
       pR[face],
       pCv[face],
       pCv[face],
       pSf[face],
       pMagSf[face]
      );
    }
    //
  }    
}





void RiemannSolver::evaluateFluxInternal
(
 scalar& rhoFlux,
 vector& rhoUFlux,
 scalar& rhoEFlux,
 const scalar& pLeft,
 const scalar& pRight,
 const vector& ULeft,
 const vector& URight,
 const scalar& TLeft,
 const scalar& TRight,
 const scalar& RLeft,
 const scalar& RRight,
 const scalar& CvLeft,
 const scalar& CvRight,
 const vector& Sf,
 const scalar& magSf
 )
{
  // Step 0 : Decode variables
  scalar rhoLeft           = pLeft/(RLeft * TLeft);
  scalar rhoRight          = pRight/(RRight * TRight);
  scalar eLeft             = CvLeft * TLeft + 0.5*magSqr(ULeft);
  scalar eRight            = CvRight * TRight + 0.5*magSqr(URight);
  const scalar kappaLeft   = (CvLeft+RLeft)/CvLeft;
  const scalar kappaRight  = (CvRight+RRight)/CvRight;
  vector normalVector      = Sf/magSf;
  const scalar contrVLeft  = (ULeft & normalVector);
  const scalar contrVRight = (URight & normalVector);
  const scalar hLeft       = eLeft + pLeft/rhoLeft;
  const scalar hRight      = eRight + pRight/rhoRight;

  // Step 1a : Compute Roe's averaged variables
  const scalar rhoTilde = sqrt(max( rhoLeft*rhoRight , SMALL ));
  const scalar rhoLeftSqrt = sqrt(max(rhoLeft,SMALL));
  const scalar rhoRightSqrt = sqrt(max(rhoRight,SMALL));
  const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
  const scalar wRight = rhoRightSqrt/(rhoLeftSqrt + rhoRightSqrt);
  const vector UTilde = ULeft*wLeft + URight*wRight;
  const scalar hTilde = hLeft*wLeft + hRight*wRight;
  const scalar qTildeSquare = magSqr(UTilde);
  const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;
  const scalar aTilde =
    sqrt(max( (kappaTilde-1)*(hTilde-0.5*qTildeSquare) ,SMALL));
  const scalar contrVTilde (UTilde & normalVector);

  // Step 1b : Compute Primitive differences
  const scalar deltaP      = pRight - pLeft;
  const scalar deltaRho    = rhoRight - rhoLeft;
  const vector deltaU      = URight - ULeft;
  const scalar deltaContrV = (deltaU & normalVector);

  // Step 4 : Compute wave strengths
  const scalar alpha1   = 
    (deltaP - rhoTilde*aTilde*deltaContrV)/(2.0*sqr(aTilde));
  const scalar alpha23 = 
    rhoTilde * ((deltaU & vector(1,1,1))  - deltaContrV);
  const scalar alpha4   = deltaRho - deltaP/sqr(aTilde);
  const scalar alpha5   =
    (deltaP + rhoTilde*aTilde*deltaContrV)/(2.0*sqr(aTilde));

  // Step 2 : Compute eigenvalues
  scalar lambda1   = mag(contrVTilde - aTilde);
  scalar lambda234 = mag(contrVTilde);
  scalar lambda5   = mag(contrVTilde + aTilde);
  //scalar lambdaMax = max(max(lambda1,lambda234),lambda5);

  // Step 3 : Compute eigenvectors
  const scalar K1_1   = 1;
  const vector K1_234 = UTilde - aTilde*normalVector;
  const scalar K1_5   = hTilde - aTilde*contrVTilde;

  const scalar K23_1   = 0;
  const vector K23_234 = vector(1,1,1) - normalVector;
  const scalar K23_5   = ( UTilde & vector(1,1,1) ) - contrVTilde ; 

  const scalar K4_1   = 1;
  const vector K4_234 = UTilde;
  const scalar K4_5   = 0.5 * qTildeSquare;

  const scalar K5_1   = 1;
  const vector K5_234 = UTilde + aTilde*normalVector;
  const scalar K5_5   = hTilde + aTilde*contrVTilde;    

  // Step 3b : Compute some alpha * eigenvector
  // For computational purpose, we will compute directly
  // alpha2*K2 + alpha3*K3
  // this is quite hardcoded, will see for a better vesrion
  const scalar alphaK23_1 = 0;
  const vector alphaK23_234 = rhoTilde * ( deltaU - deltaContrV*normalVector  );
  const scalar alphaK23_5 = rhoTilde * (( UTilde & deltaU ) - contrVTilde*deltaContrV);
  
  // Step 5 : Build actual flux
  // Step 5a : Compute Delta F
  const scalar diffF1_1   = lambda1 * alpha1 * K1_1;
  const vector diffF1_234 = lambda1 * alpha1 * K1_234;
  const scalar diffF1_5   = lambda1 * alpha1 * K1_5;
  //
  const scalar diffF23_1   = lambda234 * alphaK23_1;//alpha23 * K23_1;
  const vector diffF23_234 = lambda234 * alphaK23_234;
    //  tensor(0, alpha23,alpha23,
    //	   alpha23,0,alpha23,alpha23,alpha23,0) & K23_234;
  const scalar diffF23_5   = lambda234 * alphaK23_5;//alpha23 * K23_5;
  //
  const scalar diffF4_1   = lambda234 * alpha4 * K4_1;
  const vector diffF4_234 = lambda234 * alpha4 * K4_234;
  const scalar diffF4_5   = lambda234 * alpha4 * K4_5;
  //
  const scalar diffF5_1   = lambda5 * alpha5 * K5_1;
  const vector diffF5_234 = lambda5 * alpha5 * K5_234;
  const scalar diffF5_5   = lambda5 * alpha5 * K5_5;

  // Step 5b : Compute Left and right flux
  const scalar fluxLeft_1   = rhoLeft*contrVLeft;
  const vector fluxLeft_234 = ULeft*fluxLeft_1 + pLeft*normalVector;
  const scalar fluxLeft_5   = hLeft*fluxLeft_1;
  //
  const scalar fluxRight_1   = rhoRight*contrVRight;
  const vector fluxRight_234 = URight*fluxRight_1 + pRight*normalVector;
  const scalar fluxRight_5   = hRight*fluxRight_1;

  // Step 5c : Assemble fluxes  (Toro : (11.29))
  const scalar flux1 =
    0.5 * (fluxLeft_1 + fluxRight_1 - diffF1_1 - diffF23_1 - diffF4_1 - diffF5_1);
  const vector flux234 =
    0.5 * (fluxLeft_234 + fluxRight_234 - diffF1_234 - diffF23_234 - diffF4_234 - diffF5_234);
  const scalar flux5 =
    0.5 * (fluxLeft_5 + fluxRight_5 - diffF1_5 - diffF23_5 - diffF4_5 - diffF5_5);

  // Step 5d : Update fluxes
  rhoFlux  = flux1 * magSf;
  rhoUFlux = flux234 * magSf;
  rhoEFlux = flux5 * magSf;
}


void RiemannSolver::evaluateFluxBoundary
(
 scalar& rhoFlux,
 vector& rhoUFlux,
 scalar& rhoEFlux,
 const scalar& pLeft,
 const scalar& pRight,
 const vector& ULeft,
 const vector& URight,
 const scalar& TLeft,
 const scalar& TRight,
 const scalar& RLeft,
 const scalar& RRight,
 const scalar& CvLeft,
 const scalar& CvRight,
 const vector& Sf,
 const scalar& magSf
 )
{
  // At the moment, internal and boundary faces are evaluated
  // with the same routine ( like initial algorithm ),
  // See Jibran's for actual parsing of that..
  evaluateFluxInternal
  (
   rhoFlux,
   rhoUFlux,
   rhoEFlux,
   pLeft,
   pRight,
   ULeft,
   URight,
   TLeft,
   TRight,
   RLeft,
   RRight,
   CvLeft,
   CvRight,
   Sf,
   magSf
  );
}









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
