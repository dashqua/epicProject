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
#include "fluidThermo.H"

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
	    const scalar aTilde = Foam::sqrt(max((kappaTilde-a)*(hTilde-0.5*qTildeSquare), SMALL));
	    const scalar contrVTilde = UTilde & normalVector;

	    // Step 1b : Compute primitive differences
	    const scalar deltaP = pRight - pLeft;
	    const scalar deltaRho = rhoRight - rhoLeft;
	    const vector deltaU = URight - ULeft;

	    // Step 2 : Compute right eigenvectors
	    // note : eigen 2 and 3 are added
	    const vector l1U = UTilde - aTilde*normalVector; //tmp
	    const vector l5U = UTilde + aTilde*normalVector; //tmp  
	    const vector eigen1 = vector(1, l1U[0], l1U[1], l1U[2], hTilde-aTilde*contrVTilde);
	    const vector eigen23 = vector(0, 1, 1, 0, (UTilde & deltaU) - contrVTilde*deltaContrV);
	    const vector eigen4 = vector(1, U[0], U[1], U[2], 0.5*qTildeSquare);
	    const vector eigen5 = vector(1, l5U[0], l5U[1], l5U[2], hTilde+aTilde*contrVTilde);

	    // Step 3 : Compute eigenvalues
	    scalar lambda1 = mag(contrVTilde - aTilde);
	    scalar lambda2 = mag(contrVTilde);
	    scalar lambda3 = mag(contrVTilde + aTilde);
	    scalar lambdaMax = max(max(lambda1,lambda2),lambda3);

	    // Step 4 : Compute wave strengths
	    const scalar r1 = (deltaP - rhoTilde*deltaContrV)/(2*Foam::sqrt(aTilde));
	    //const scalar r23 = /!\ HAS TO BE CODED
	    const scalar r4 = deltaRho - deltaP/Foam::sqrt(aTilde);
	    const scalar r5 = (deltaP + rhoTilde*deltaContrV)/(2*Foam::sqrt(aTilde));

	    // Step 5 : Assemble the flux
	    // Step 5a : Compute flux differences
	    const vector diffF1 = lambdaMax*r1*eigen1;
	    const vector diffF23 = lambdaMax*(r23*eigen23);
	    const vector diffF4 = lambdaMax*r4*eigen4;
	    const vector diffF5 = lambdaMax*r5*eigen5;
	    // Step 5b : Compute left and right fluxes
	    const vector urcvL = ULeft*rhoLeft*contrVLeft;// tmp
	    const vector fluxLeft1 = vector(rhoLeft*contrVLeft,
					    urcvL[0],
					    urcv[1],
					    urcv[2],
					    hLeft*rhoLeft*contrVLeft
					    );  // REECRIRE CE VECTOR + FAIRE LES AUTRES
	  }
	
	runTime.write();
      
	Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;

      }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
