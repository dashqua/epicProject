/*---------------------------------------------------------------------------* \
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

Application
    incredibleSolver

Description
    Density-based compressible explicit time-marching flow solver
    using enthalpy based thermo packages

Author
    Thomas Di Giusto, based on the work of Hrvoje Jasak and Aleksandar Jemcov

\*---------------------------------------------------------------------------*/

// Core headers
#include "fvCFD.H"
#include "psiThermo.H"
#include "bound.H"

// Limiters
#include "MDLimiter.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"

// Custom libraries
#include "arbMesh.H"
#include "RiemannSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{  
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createTimeControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    Info<< "\nStarting time loop\n" << endl;

    // Runge-Kutta coefficient
    scalarList beta(4);
    beta[0] = 0.1100;
    beta[1] = 0.2766;
    beta[2] = 0.5000;
    beta[3] = 1.0000;

    // Switch off solver messages
    lduMatrix::debug = 0;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "readFieldBounds.H"

//////// The following replaces #include "compressibleCourantNo.H"
      scalar CoNum = 0.0;
      scalar meanCoNum = 0.0;
      if (mesh.nInternalFaces())
      {
	volScalarField speed_of_sound = sqrt(thermo->Cp() / thermo->Cv() * (thermo->Cp() - thermo->Cv()) * T);
	surfaceScalarField acCo =
	  (mag(phi) / (fvc::interpolate(rho) * mesh.magSf()) + fvc::interpolate(speed_of_sound))
	  * mesh.surfaceInterpolation::deltaCoeffs() * runTime.deltaT();
	
	CoNum = gMax(acCo.internalField());

	meanCoNum = gSum(acCo.internalField())  / mesh.nInternalFaces();
      }
      Info << "Acoustic Courant Number mean: " << meanCoNum
	   << " max: " << CoNum << endl;
//////// 
#       include "setDeltaT.H"

        runTime++;

        Info<< "\n Time = " << runTime.value() << endl;

        // Low storage Runge-Kutta time integration
        forAll (beta, i)
        {
	    // Compute EUL variables from TALE variables
	    //aMsh.computeEULfromTALE();
	  
	    // Solve the approximate Riemann problem for this time step
	    RS.computeFlux( rhoFlux, rhoUFlux, rhoEFlux, gradP, gradU, gradT);

	    // Use EUL variables to get TALE variables
	    //aMsh.computeTALEfromEUL();

            // Time integration
            solve
            (
                1.0/beta[i]*fvm::ddt(rho)
              + fvc::div(rhoFlux)
            );

            solve
            (
                1.0/beta[i]*fvm::ddt(rhoU)
              + fvc::div(rhoUFlux)
            );

            solve
            (
                1.0/beta[i]*fvm::ddt(rhoE)
              + fvc::div(rhoEFlux)
            );

#           include "updateFields.H"

        }

	// theoretical variables are updated.
	//aMsh.updateFields();	
	
        runTime.write();

        Info<< "    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl << endl << endl << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
