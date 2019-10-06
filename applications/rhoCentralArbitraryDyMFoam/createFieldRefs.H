//volScalarField& p = thermo.p();
volScalarField p 
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//thermo.rho()                                                                                                                                                                                                                        
);



//const volScalarField& T = thermo.T();
const volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh//thermo.rho()                                                                                                                                                                                                                        
);

const volScalarField psi
(
 IOobject
 (
  "psi",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("psi",dimensionSet(0, -2, 2, 0, 0), scalar(1./(0.08205*300.)))
 );



//const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}
