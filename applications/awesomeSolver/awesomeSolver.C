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
    awesomeSolver

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop \n" << endl;

    // Choose Roe's original eigenvectors (1) or its light version (2)
    const bool Roe_eigenvectors = 0;    
    
    while (runTime.run())
      {
	Info << "Time = " << runTime.timeName() << nl << endl;

	// CODE //
	//COMPUTE THE FLUX
	const unallocLabelList& owner = mesh.owner();
	const unallocLabelList& neighbour = mesh.neighbour();
	const surfaceVectorField& Sf = mesh.Sf();
	const surfaceScalarField& magSf = mesh.magSf();
	const volScalarField Cv = thermo.Cv();
	const volScalarField R = thermo.Cp() - Cv;
	gradP = fvc::grad(p);
	//gradP.correctBoundaryConditions;
	gradU = fvc::grad(U);
	//gradU.correctBoundaryConditions;
	const volVectorField& cellCentre =mesh.C();
	const surfaceVectorField& faceCentre = mesh.Cf();
	
	forAll(owner, faceI)
	  {
	    const label own = owner[faceI];
	    const label nei = neighbour[faceI];
	    const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
	    const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

	    //==>
	    // Step 0 : Decode left and right info 
	    vector normalVector = Sf[faceI]/magSf[faceI];
	    scalar rhoLeft  = rho[own];
	    scalar rhoRight = rho[nei];
	    scalar eLeft  = E[own];
	    scalar eRight = E[nei];
	    const scalar pLeft  = p[own];
	    const scalar pRight = p[nei];
	    const vector URight = U[nei];
	    const vector ULeft  = U[own];
	    const scalar kappaLeft  = (Cv[own]+R[own])/Cv[own];
	    const scalar kappaRight = (Cv[nei]+R[nei])/Cv[nei];
	    const scalar contrVLeft  = U[own] & normalVector;
	    const scalar contrVRight = U[nei] & normalVector;
	    const scalar hLeft  = eLeft + pLeft/rhoLeft;
	    const scalar hRight = eRight + pRight/rhoRight;

	    // Step 1a : Compute Roe's averaged quantities
	    const scalar rhoTilde = Foam::sqrt(max(rhoLeft*rhoRight, SMALL));
	    const scalar rhoLeftSqrt  = Foam::sqrt(max(rhoLeft, SMALL));
	    const scalar rhoRightSqrt = Foam::sqrt(max(rhoRight, SMALL));
	    const scalar wLeft  = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
	    const scalar wRight = 1 - wLeft;
	    const vector UTilde = ULeft*wLeft + URight*wRight;
	    const scalar hTilde = hLeft*wLeft + hRight*wRight;
	    const scalar qTildeSquare = magSqr(UTilde);
	    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;
	    const scalar aTilde = Foam::sqrt(max((kappaTilde-1)*(hTilde-0.5*qTildeSquare), SMALL));
	    const scalar contrVTilde = UTilde & normalVector;

	    // Step 1b : Compute primitive differences
	    const scalar deltaP      = pRight - pLeft;
	    const scalar deltaRho    = rhoRight - rhoLeft;
	    const vector deltaU      = URight - ULeft;
	    const scalar deltaContrV = (deltaU & normalVector);
	    
	    // Step 2 : Compute right eigenvectors
	    // note : eigen 2 and 3 are added
	    // note :  5-components vectors have to be separated because OF does not have that structure
	    // so eigen_i_j is the ij j-th coef of the i-th eigenvector, or the ij coef of the P matrix
	    // note : for convenience, eigen_i_234 refers to 2,3, and 4th coeffs of the i-th eigenvector    
	    // first eigenvect is vector(1, u-a, v, w, H-ua);
	    const scalar eigen1_1   = 1;
	    const vector eigen1_234 = UTilde - aTilde*normalVector;
	    const scalar eigen1_5   = hTilde-aTilde*contrVTilde;
	    if (Roe_eigenvectors) // Roe eigenvectors
	      {  // note : NOT CODED YET (see 5th component)
	    // second+third eigenvect is vector(0, 0, v, w, v^2 + w^2);
	    const scalar eigen23_1   = 0;
	    const vector eigen23_234 = deltaU - deltaContrV*normalVector;
	    const scalar eigen23_5   = (UTilde & deltaU) - contrVTilde*deltaContrV;
	      }
	    else // lightRoe eigenvectors
	      {
	    // second+third eigenvect is vector(0, 0, 1, 1, v + w);
	    const scalar eigen23_1   = 0;
	    const vector eigen23_234 = vector(1,1,1) - normalVector;
	    const scalar eigen23_5   = (UTilde & vector(1,1,1)) - contrVTilde;
	      }
	    const vector eigen23_234bis = vector(1,1,1) - normalVector;
	    // fourth eigenvect is vector(1, u, v, w, 0.5*q^2)
	    const scalar eigen4_1   = 1;
	    const vector eigen4_234 = UTilde;
	    const scalar eigen4_5   = .5 * qTildeSquare;
	    // fifth eigenvect is vector(1, u+a, v, w, H+au)
	    const scalar eigen5_1   = 1;
	    const vector eigen5_234 = UTilde + aTilde*normalVector;
	    const scalar eigen5_5   = hTilde-aTilde*contrVTilde;
	    //const vector eigen5 = vector(1, l5U[0], l5U[1], l5U[2], hTilde+aTilde*contrVTilde);

	    // Step 3 : Compute eigenvalues
	    scalar lambda1 = mag(contrVTilde - aTilde);
	    scalar lambda234 = mag(contrVTilde);
	    scalar lambda5 = mag(contrVTilde + aTilde);
	    scalar lambdaMax = max(max(lambda1,lambda234),lambda5);

	    // Step 4 : Compute wave strengths
	    const scalar a1  = (deltaP - rhoTilde*deltaContrV)/(2*Foam::sqrt(aTilde));
	    const scalar a23 = rhoTilde*((deltaU - deltaContrV*normalVector) & vector(1,1,1));
	    const scalar a4  = deltaRho - deltaP/Foam::sqrt(aTilde);
	    const scalar a5  = (deltaP + rhoTilde*deltaContrV)/(2*Foam::sqrt(aTilde));

	    // Step 5 : Assemble the flux
	    // Step 5a : Compute left and right fluxes
	    const scalar upvarphiLeft  = rhoLeft * contrVLeft; // tmp
	    const scalar upvarphiRight = rhoRight * contrVRight; // tmp
	    // note : like eigenvectors, the fluxes are 5-vectors	    
	    const scalar fluxLeft1   = upvarphiLeft;
	    const vector fluxLeft234 = ULeft * upvarphiLeft + pLeft*normalVector;
	    const scalar fluxLeft5   = upvarphiLeft * hLeft;
	    //
	    const scalar fluxRight1   = upvarphiRight;
	    const vector fluxRight234 = URight * upvarphiRight + pRight*normalVector;
	    const scalar fluxRight5   = upvarphiRight * hRight;

	    // Step 5b : Compute flux differences or deltas
	    // note : Compute the 5 flux deltas
	    // They are then added to be the stabilization term
	    const scalar diffF1_1   = lambda1 * a1 * eigen1_1;
	    const vector diffF1_234 = lambda1 * a1 * eigen1_234;
	    const scalar diffF1_5   = lambda1 * a1 * eigen1_5;
	    //
	    // note : see the doc for explanation on waveTensor
	    // note : can obviously be improved
	    const scalar a2 = a23, a3 = a23;
	    tensor waveTensor( 0, a2, a3,
			      a2,  0, a3,
			      a2, a3,  0);
	    vector waveTensorProjection = waveTensor & normalVector; // (0,a2,a3) for x-sweep
	    const scalar diffF23_1   = 0; // derived by hand	    
	    const vector diffF23_234 = lambda234 * //waveTensor * eigen23_234;
 	      vector(waveTensorProjection[0] * eigen23_234bis[0],
		     waveTensorProjection[1] * eigen23_234bis[1],    // BIS A ENLEVER
		     waveTensorProjection[2] * eigen23_234bis[2]);
	    const scalar diffF23_5   = lambda234 *
	      vector(waveTensorProjection[0] * UTilde[0] * deltaU[0],
		     waveTensorProjection[1] * UTilde[1] * deltaU[1],
		     waveTensorProjection[2] * UTilde[2] * deltaU[2])
	      & vector(1,1,1);

	    const scalar diffF4_1   = lambda234 * a4 * eigen4_1;
	    const vector diffF4_234 = lambda234 * a4 * eigen4_234;
	    const scalar diffF4_5   = lambda234 * a4 * eigen4_5;

	    const scalar diffF5_1   = lambda5 * a5 * eigen5_1;
	    const vector diffF5_234 = lambda5 * a5 * eigen5_234;
	    const scalar diffF5_5   = lambda5 * a5 * eigen5_5;

	    // step 5b : Assembly of the physical fluxes
	    const scalar flux_1    = 0.5 *
	      (fluxLeft1 + fluxRight1 - diffF1_1 - diffF23_1 - diffF4_1 - diffF5_1); 
	    const vector flux_234  = 0.5 *
	      (fluxLeft234 + fluxRight234 - diffF1_234 - diffF23_234 - diffF4_234 - diffF5_234);
	    const scalar flux_5    = 0.5 *
	      (fluxLeft5 + fluxRight5 - diffF1_5 - diffF23_5 - diffF4_5 - diffF5_5);


	    rhoFlux[faceI]  = flux_1   * magSf[faceI];
	    rhoUFlux[faceI] = flux_234 * magSf[faceI];
	    EFlux[faceI]    = flux_5   * magSf[faceI];

	      // note : most of the calculation can be vectorized
	  }

	// Time integration
	solve(fvm::ddt(rho)  + fvc::div(rhoFlux));
	solve(fvm::ddt(rhoU) + fvc::div(rhoUFlux));
	solve(fvm::ddt(E)    + fvc::div(EFlux));
	
	// Update Fields
	U = rhoU / U;
	// U.correctBoundaryConditions();
	H = Cp/Cv*(E - 0.5*magSqr(U));
	//H.correctBoundaryConditions();
	dimensionedScalar CpMin = min(Cp);
	dimensionedScalar CpMax = max(Cp);
	dimensionedScalar hMin = CpMin*TMin;
	dimensionedScalar hMax = CpMax*TMax;








	
	runTime.write();
      
	Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;

      }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
