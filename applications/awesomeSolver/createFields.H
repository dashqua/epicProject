Info << "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
 fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");  // probablement enables  h and e fields to be taken from thermo

Info << "Reading pressure p and enthalpy h\n" << endl;
volScalarField& p = thermo.p();
volScalarField& h = thermo.he();

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
 p/(thermo.gamma()-1) + 0.5 * magSqr(U)
);

Info << "Reading Flux Fields\n" << endl;
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
 linearInterpolate(p/(thermo.gamma()-1 + 0.5*magSqr(U)))
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
