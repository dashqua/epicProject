/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PRE_rusanovFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PRE_rusanovFlux::evaluateFlux
(
    scalar& rhoFlux_TALE,
    vector& rhoUFlux_TALE,
    scalar& rhoEFlux_TALE,
    const scalar& pLeft_TALE,
    const scalar& pRight_TALE,
    const vector& ULeft_TALE,
    const vector& URight_TALE,
    const scalar& TLeft_TALE,
    const scalar& TRight_TALE,
    const scalar& RLeft_TALE,
    const scalar& RRight_TALE,
    const scalar& CvLeft_TALE,
    const scalar& CvRight_TALE,
    const vector& Sf,
    const scalar& magSf,
    vector& xyzOwn,
    vector& xyzNei,
    const double& t,
    const scalar& magSfOld,
    const objectRegistry& db,
    arbMesh& amsh,
    label& faceI
) const
{
  // Step 0: Conversion from U_TALE to U_E
  // Note : Own <=> Left ; Nei <=> Right
  // /!\ Flux are interpolated before in numerifFlux
  scalar xLeft = xyzOwn[0];
  scalar yLeft = xyzOwn[1];
  scalar xRight = xyzNei[0];
  scalar yRight = xyzNei[1];
  scalar jwLeft = amsh.jw(xyzOwn);
  scalar jwRight = amsh.jw(xyzNei);
  
  //scalar& rhoFlux = rhoFlux_TALE / jw;
  //vector& rhoUFlux_TALE;
  //scalar& rhoEFlux_TALE;
  const scalar& pLeft = pLeft_TALE / jwLeft;
  const scalar& pRight = pRight_TALE / jwRight;
  const vector& ULeft = ULeft_TALE / jwLeft;
  const vector& URight = URight_TALE / jwRight;
  const scalar& TLeft = TLeft_TALE/ jwLeft;
  const scalar& TRight = TRight_TALE / jwRight;
  const scalar& RLeft = RLeft_TALE / jwLeft;
  const scalar& RRight = RRight_TALE / jwRight;
  const scalar& CvLeft = CvLeft_TALE / jwLeft;
  const scalar& CvRight = CvRight_TALE / jwRight;
  
  // Step 1: decode rho left and right:
    scalar rhoLeft = pLeft/(RLeft*TLeft);
    scalar rhoRight = pRight/(RRight*TRight);

    // Decode left and right total energy:
    scalar eLeft = CvLeft*TLeft+0.5*magSqr(ULeft);
    scalar eRight = CvRight*TRight+0.5*magSqr(URight);

    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
    const scalar kappaLeft = (CvLeft+RLeft)/CvLeft;
    const scalar kappaRight = (CvRight+RRight)/CvRight;

    // normal vector
    vector normalVector = Sf/magSf;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Compute left and right total enthalpies:
    const scalar hLeft = eLeft + pLeft/rhoLeft;
    const scalar hRight = eRight + pRight/rhoRight;

    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = hLeft*wLeft + hRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Speed of sound
    const scalar cTilde =
        sqrt(max((kappaTilde - 1)*(hTilde - 0.5*qTildeSquare), SMALL));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

    // Roe and Pike - formulation
    const scalar r1 =
        (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
    const scalar r2 = deltaRho - deltaP/sqr(cTilde);
    const scalar r3 =
        (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // Step 5: compute l vectors

    // rho row:
    const scalar l1rho = 1;
    const scalar l2rho = 1;
    const scalar l3rho = 0;
    const scalar l4rho = 1;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;

    // second U column
    const vector l2U = UTilde;

    // third U column
    const vector l3U = deltaU - deltaContrV*normalVector;

    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // Step 6a: compute eigenvalues

    // derived from algebra by hand, only for Euler equation usefull
    scalar lambda1 = mag(contrVTilde - cTilde);
    scalar lambda2 = mag(contrVTilde);
    scalar lambda3 = mag(contrVTilde + cTilde);
    
    // Step 6b: Shift with mapping coefficient
    //#   include "createShiftFields.H"
    //#   include "mappingShift.H"
    amsh.Shift(lambda1, xyzOwn, xyzNei, faceI, Sf/magSf);
    amsh.Shift(lambda1, xyzOwn, xyzNei, faceI, Sf/magSf);
    amsh.Shift(lambda1, xyzOwn, xyzNei, faceI, Sf/magSf);
    
    scalar lambdaMax = max(max(lambda1,lambda2),lambda3);

    // Step 7: Compute flux differences

    // Components of deltaF1
    const scalar diffF11 = lambdaMax*r1*l1rho;
    const vector diffF124 = lambdaMax*r1*l1U;
    const scalar diffF15 = lambdaMax*r1*l1e;

    // Components of deltaF2
    const scalar diffF21 = lambdaMax*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambdaMax*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambdaMax*(r2*l2e + rhoTilde*l3e);

    // Components of deltaF3
    const scalar diffF31 = lambdaMax*r3*l4rho;
    const vector diffF324 = lambdaMax*r3*l4U;
    const scalar diffF35 = lambdaMax*r3*l4e;

    // Step 8: compute left and right fluxes

    // Left flux 5-vector
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = hLeft*fluxLeft11;

    // Right flux 5-vector
    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = hRight*fluxRight11;

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxLeft11 + fluxRight11 - (diffF11 + diffF21 + diffF31));
    const vector flux24 =
        0.5*(fluxLeft124 + fluxRight124 - (diffF124 + diffF224 + diffF324));
    const scalar flux5 =
        0.5*(fluxLeft15 + fluxRight15 - (diffF15 + diffF25 + diffF35));

    // Compute private data
    // NOTE : this is F_E
    scalar rhoFlux_E  = flux1*magSf;
    vector rhoUFlux_E = flux24*magSf;
    scalar rhoEFlux_E = flux5*magSf;

    // Conversion back to F_TALE
    //rhoFlux = (amsh.deltaw(xRight,yRight)+amsh.deltaw(xLeft,yLeft))/2 * ( rhoFlux_E - (amsh.Uwn(xRight,yRight)+amsh.Uwn(xLeft,yLeft))/2 *   U );
    //rhoUFlux = rhoUFlux_E;
    //rhoEFlux = RhoEFlux_E;
    
    // the physical flux is sum of the left and right numerical fluxes
    // in addition to a correction term (stabilization) 
    rhoFlux_TALE = \
      ( amsh.deltaw(xyzNei, faceI) * (rhoFlux_E - amsh.Uwn(xyzNei, Sf/magSf) * rhoRight) + \
        amsh.deltaw(xyzOwn, faceI) * (rhoFlux_E - amsh.Uwn(xyzOwn, Sf/magSf) * rhoRight) ) /2;
    //
    rhoUFlux_TALE = \
      ( amsh.deltaw(xyzNei, faceI) * (rhoUFlux_E - amsh.Uwn(xyzNei, Sf/magSf) * URight) + \
        amsh.deltaw(xyzOwn, faceI) * (rhoUFlux_E - amsh.Uwn(xyzOwn, Sf/magSf) * ULeft) ) /2;
    //
    rhoEFlux_TALE = \
      ( amsh.deltaw(xyzNei, faceI) * (rhoEFlux_E - amsh.Uwn(xyzNei, Sf/magSf) * rhoRight*eRight) + \
        amsh.deltaw(xyzOwn, faceI) * (rhoEFlux_E - amsh.Uwn(xyzOwn, Sf/magSf) * rhoLeft*eLeft) ) /2;
    //stabilization here;;
   
}

// ************************************************************************* //
