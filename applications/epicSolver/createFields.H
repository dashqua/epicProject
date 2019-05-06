    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> thermo
    (
        basicPsiThermo::New(mesh)
    );

    // Primitive variables

    volScalarField& h = thermo->h();
    volScalarField& p = thermo->p();
    const volScalarField& T = thermo->T();

    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo->rho()
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Material (cell centered) x
    const volVectorField& C = mesh.C();
    volVectorField x
    (
        IOobject("x", mesh),
	C
     );

    // Nodal displacements
    pointMesh pMesh(mesh);
    pointVectorField xN
    (
        IOobject("xN", mesh),
	pMesh,
	dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    );
    xN.internalField() = mesh.points();
    pointVectorField XN = xN;
    pointVectorField uN
    (
        IOobject
        (
	    "uN",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	pMesh,
	dimensionedVector("uN", dimLength, vector::zero)
    );
    uN.write();

    // Normals
    const surfaceVectorField& Sf = mesh.Sf();
    //const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField N = Sf/mesh.magSf();
    surfaceVectorField n
    (
        IOobject("n", mesh),
	N
    );

    // Conservative variables

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    volPointInterpolation interpolateVolPoint(mesh);
    pointVectorField UN = interpolateVolPoint.interpolate(U);
/*
    pointVectorField rhoUN
    (
        IOobject
        (
            "rhoUN",
            runTime.timeName(),
            pMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	rho*UN
    );
*/
/*
    pointVectorField UN
    (
        IOobject
        (
            "UN",
            runTime.timeName(),
            pMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh
    );
*/

    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(h + 0.5*magSqr(U)) - p
    );


    // Create numeric flux
    autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New
    (
        p,
        U,
        T,
        thermo()
    );
    basicNumericFlux& dbnsFlux = dbnsFluxPtr();

    // Create mass flux alias for easier coupling with other code components
    const surfaceScalarField& phi = dbnsFlux.rhoFlux();

    // delta T
    dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());

    // interpolation schemes
    //interpolationSchemes interpolate(mesh);
//fvc::interpolate.volToPoint(rhoU, rhoUGrad, rhoUN)

    // mesh displacement vector
    pointVectorField MDN
    (
        IOobject
        (
	    "MDN",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	pMesh,
	dimensionedVector("MDN", dimLength, vector::zero)
    );
    MDN.write();
