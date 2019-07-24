diff --git a/createFields.H b/home/thomas/foam/foam-extend-4.0/applications/solvers/compressible/dbnsFoam/createFields.H
index b6a9974..8b0f268 100644
--- a/createFields.H
+++ b/home/thomas/foam/foam-extend-4.0/applications/solvers/compressible/dbnsFoam/createFields.H
@@ -33,61 +33,12 @@
             "U",
             runTime.timeName(),
             mesh,
-            IOobject::NO_READ, // MUST_READ,
+            IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
-        mesh,
-	dimensionedVector
-	(
-	    "U",         
-	    dimensionSet(0,1,-1,0,0,0,0),// Foam::sin(PI*C.component(1)/(2*height));
-	    Foam::vector(1 * Foam::cos(26.56), 1 * Foam::sin(26.56), 0)
-	)
-     );        //  	dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
-
-
-    // Material (cell centered) x
-    const volVectorField& C = mesh.C();
-    volVectorField x
-    (
-        IOobject("x", mesh),
-	C
-     );
-
-    // Nodal displacements
-    pointMesh pMesh(mesh);
-    pointVectorField xN
-    (
-        IOobject("xN", mesh),
-	pMesh,
-	dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
-    );
-    xN.internalField() = mesh.points();
-    pointVectorField XN = xN;
-    pointVectorField uN
-    (
-        IOobject
-        (
-	    "uN",
-	    runTime.timeName(),
-	    mesh,
-	    IOobject::NO_READ,
-	    IOobject::AUTO_WRITE
-	),
-	pMesh,
-	dimensionedVector("uN", dimLength, vector::zero)
+        mesh
     );
-    uN.write();
 
-    // Normals
-    const surfaceVectorField& Sf = mesh.Sf();
-    //const surfaceScalarField& magSf = mesh.magSf();
-    const surfaceVectorField N = Sf/mesh.magSf();
-    surfaceVectorField n
-    (
-        IOobject("n", mesh),
-	N
-    );
 
     // Conservative variables
 
@@ -104,37 +55,6 @@
         rho*U
     );
 
-    volPointInterpolation interpolateVolPoint(mesh);
-    pointVectorField UN = interpolateVolPoint.interpolate(U);
-/*
-    pointVectorField rhoUN
-    (
-        IOobject
-        (
-            "rhoUN",
-            runTime.timeName(),
-            pMesh,
-            IOobject::NO_READ,
-            IOobject::NO_WRITE
-        ),
-	rho*UN
-    );
-*/
-/*
-    pointVectorField UN
-    (
-        IOobject
-        (
-            "UN",
-            runTime.timeName(),
-            pMesh,
-            IOobject::NO_READ,
-            IOobject::NO_WRITE
-        ),
-        pMesh
-    );
-*/
-
     volScalarField rhoE
     (
         IOobject
@@ -161,86 +81,3 @@
 
     // Create mass flux alias for easier coupling with other code components
     const surfaceScalarField& phi = dbnsFlux.rhoFlux();
-
-    // delta T
-    dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());
-
-    // interpolation schemes
-    //interpolationSchemes interpolate(mesh);
-//fvc::interpolate.volToPoint(rhoU, rhoUGrad, rhoUN)
-
-    // mesh displacement vector
-    Info << "Reading field MDN\n" << endl;
-    pointVectorField MDN
-    (
-        IOobject
-        (
-	    "MDN",
-	    runTime.timeName(),
-	    mesh,
-	    IOobject::NO_READ,
-	    IOobject::AUTO_WRITE
-	),
-	pMesh,
-	dimensionedVector("MDN", dimLength, vector::zero)
-    );
-    MDN.write();
-
-    volScalarField JW
-    (
-     IOobject
-     ("JW_",
-      runTime.timeName(),
-      mesh,
-      IOobject::NO_READ,
-      IOobject::AUTO_WRITE
-      ),
-     (1 * rho) - thermo->rho()
-     );
-
-    Info<< "Reading analytic fields\n" << endl;
-    volVectorField U_theo
-    (
-        IOobject
-        (
-            "U",
-            runTime.timeName(),
-            mesh,
-            IOobject::NO_READ, // MUST_READ,
-            IOobject::AUTO_WRITE
-        ),
-        mesh,
-	dimensionedVector
-	(
-	    "U_theo",         
-	    dimensionSet(0,1,-1,0,0,0,0),// Foam::sin(PI*C.component(1)/(2*height));
-	    Foam::vector(1 * Foam::cos(26.56), 1 * Foam::sin(26.56), 0)
-	)
-     );        //  	dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
-
-    volScalarField rho_theo
-    (
-        IOobject
-        (
-            "rho_theo",
-            runTime.timeName(),
-            mesh,
-            IOobject::NO_READ,
-            IOobject::AUTO_WRITE
-        ),
-        thermo->rho()
-    );
-
-    volScalarField p_theo
-    (
-        IOobject
-        (
-            "p_theo",
-            runTime.timeName(),
-            mesh,
-            IOobject::NO_READ,
-            IOobject::AUTO_WRITE
-        ),
-        3.0*rho/rho
-    );
-
