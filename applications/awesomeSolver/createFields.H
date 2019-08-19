Info << "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
 psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
//thermo.validate(args.executable(), "h", "e");  // probablement enables  h and e fields to be taken from thermo

Info << "Reading pressure p, internal energy e\n" << endl;
volScalarField& p = thermo.p();
volScalarField& e = thermo.he(); // TO VERIFY

Info << "Reading capcities Cp and Cv\n" << endl;
volScalarField& Cp = thermo.Cp();
volScalarField& Cv = thermo.Cv();

volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 thermo.rho()
);

Info << "Reading field U\n" << endl;
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

Info << "Reading field H\n" << endl;
volScalarField H
(
 IOobject
 (
  "H",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 p*(thermo.gamma()/(thermo.gamma()-1))/rho + 0.5*magSqr(U)
);

Info << "Reading field E\n" << endl;
volScalarField E
(
 IOobject
 (
  "E",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 rho*H - p//p/(thermo.gamma()-1) + 0.5 * rho * magSqr(U) // e
);

Info << "Reading conserved variable rhoU\n" << endl;
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


Info << "Reading rho Flux Fields\n" << endl;
surfaceScalarField rhoFlux
(
 IOobject
 (
  "phi",
  mesh.time().timeName(),mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 linearInterpolate(thermo.rho()*U) & mesh.Sf()
);

Info << "Reading rho*U Flux Fields\n" << endl;
surfaceVectorField rhoUFlux
(
 IOobject
 (
  "rhoUFlux",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 rhoFlux*linearInterpolate(U)
);

Info << "Reading E Flux Fields\n" << endl;
surfaceScalarField EFlux
(
 IOobject
 (
  "EFlux",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 linearInterpolate(U*(E+p)) & mesh.Sf() //linearInterpolate( U*(E+p) )   //linearInterpolate(p/(thermo.gamma()-1 + 0.5*magSqr(U)))
);

Info << "Reading Gradient Fields\n" << endl;
volVectorField gradP(fvc::grad(p));
volTensorField gradU(fvc::grad(U));




Info << "Creating Phi\n" << endl;
#include "compressibleCreatePhi.H"

//Info << "Initiating pressureControl\n" << endl;
//pressureControl pressureControl(p, rho, simple.dict());

// ??
mesh.setFluxRequired(p.name());

Info << "Initializing initialMassan" << endl;
dimensionedScalar initialMass = fvc::domainIntegrate(rho);


#include "createMRF.H"
//#include "createFvOptions.H"
