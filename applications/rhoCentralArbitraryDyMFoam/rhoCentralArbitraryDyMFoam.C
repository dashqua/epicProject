/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    rhoCentralDyMFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.update();

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));
        rho.correctBoundaryConditions();
	
        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() ==
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        //rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
	      - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
	     + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
	  - fvc::div(sigmaDotU)
        );

        e == rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        //thermo.correct();
	/*
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
	*/
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
            );
            //thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }

	///////////////////////////////
	/// ENFORCEMENT OF PRESSURE ///
	///////////////////////////////
	scalar M  = 0.5;
	scalar M2 = M*M;
	scalar I  = 5.0;
	scalar I2 = I*I;
	scalar r  = 1.5;
	scalar r2 = r*r;
	scalar theta = Foam::atan(0.5);
	//scalar Rhoinf = 1;
	scalar Uinf = 0.8944;
	scalar Vinf = 0.4472;
	scalar Pinf = 3.;
	scalar pii  = constant::mathematical::pi;
	scalar pii2 = pii*pii;
	scalar t = mesh.time().value();
	volVectorField C = mesh.C();	
	forAll(C, cell)
	  {
	    scalar gamma = 1.3;//Cp_[cell]/Cv_[cell];                                                                                                                        
	    scalar x1 = C[cell].x();
	    scalar x2 = C[cell].y();
	    scalar v1 = Uinf*Foam::cos(theta);
	    scalar v2 = Vinf*Foam::sin(theta);
	    p[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , gamma/(gamma-1)) ;
	  }

	p.correctBoundaryConditions();
	Info << p[2782] <<endl;
	///////////////////////////////
	/// ENFORCEMENT OF BOUNDARY ///
	///////////////////////////////
	/*	
	const polyBoundaryMesh& bm = mesh.boundaryMesh();
	const label& sidesPatchID = bm.findPatchID("emptyPatches_empt"); 
	p.boundaryFieldRef().set(sidesPatchID,fvPatchField<scalar>::New("empty", mesh.boundary()[sidesPatchID], p));
	const label& topPatchID = bm.findPatchID("top_cyc");
	p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[topPatchID], p));
	const label& bottomPatchID = bm.findPatchID("bottom_cyc");
	p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[bottomPatchID], p));
	const label& inletPatchID = bm.findPatchID("inlet_cyc");
	p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[inletPatchID], p));
	const label& outletPatchID = bm.findPatchID("outlet_cyc");
	p.boundaryFieldRef().set(topPatchID,fvPatchField<scalar>::New("fixedValue", mesh.boundary()[outletPatchID], p));
	*/
	/*
	U.boundaryFieldRef().set(sidesPatchID,fvPatchField<vector>::New("empty", mesh.boundary()[sidesPatchID], U));
	U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[topPatchID], U));
	U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[bottomPatchID], U));
	U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[inletPatchID], U));
	U.boundaryFieldRef().set(topPatchID,fvPatchField<vector>::New("fixedValue", mesh.boundary()[outletPatchID], U));	
	*/
	/*
        p.ref() =  
	  rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
	*/

	///////////////////////////////
	///// END OF ENFORCEMENT //////
	///////////////////////////////

	
        turbulence->correct();

        //runTime.write();
	if (runTime.writeTime() )
	  {
	    U.write();
	    T.write();
	    rho.write();
	    p.write();
	  }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
