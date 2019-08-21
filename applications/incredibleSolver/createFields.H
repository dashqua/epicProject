Info << "Reading thermophysical properties\n" << endl;
autoPtr<basicPsiThermo> thermo
(
 basicPsiThermo::New(mesh)
 );
// Primitive variables
Info << "\nReading primitive variables\n";
Info << " [+] Reading enthalpy h\n";
volScalarField& h = thermo->h();
Info << " [+] Reading pressure p\n";
volScalarField p
(
 IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
 mesh // thermo->p()
);
Info << " [+] Reading temperature T\n";
const volScalarField& T = thermo->T();
Info << " [+] Reading density rho\n";
volScalarField rho
(
 IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 thermo->rho()
);
Info << " [+] Reading velocity U\n";
volVectorField U
(
 IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
 mesh
);
Info << "\nReading arbitrary Mesh variables\n";
Info << " [+] Reading cell centers\n";
const volVectorField& C = mesh.C();
Info << " [+] Reading normals N and n\n";
const surfaceVectorField& Sf = mesh.Sf();
//const surfaceVectorField N = Sf/mesh.magSf();
//surfaceVectorField n(IOobject("n", mesh),N);
Info << " [+] Reading pMesh\n";
pointMesh pMesh(mesh);
//Info << " [+] Reading current position x\n";
//volVectorField x(IOobject("x",mesh), C);
Info << " [+] Reading Mesh Displacement Vector\n";
pointVectorField MDN
(
 IOobject("MDN", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 pMesh,
 dimensionedVector("MDN", dimLength, vector::zero)
);
Info << " [+] Reading arbitrary determinant";
volScalarField JW
(
 IOobject("JW", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 0 * p
);

Info << "\nReading TALE fields\n";
Info << " [+] Reading rho_TALE\n";
volScalarField rho_TALE
(
 IOobject
 (
  "rho_TALE",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 thermo->rho()
);

Info << " [+] Reading U_TALE\n";
volVectorField U_TALE
(
 IOobject
 (
  "U_TALE",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 U
);


Info << "\nReading analytic fields\n";
Info << " [+] Reading theoretical U\n";
volVectorField U_theo
(
 IOobject("U_theo", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 mesh,
 dimensionedVector("U_theo", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0, 0, 0) )
);

Info << " [+] Reading theoretical rho\n";
volScalarField rho_theo
(
 IOobject("rho_theo", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 3 * (rho/rho)
 );

Info << " [+] Reading theoretical p\n";
volScalarField p_theo
(
 IOobject("p_theo", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 3 * (rho/rho)
);


// Conservative variables
Info << "\nReading Conservative variables\n";
Info << " [+] Reading rho*U\n";
volVectorField rhoU
(
 IOobject("rhoU", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 rho*U
);
Info << " [+] Reading rho*E";
volScalarField rhoE
(
 IOobject("rhoE", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 rho*h - p // NOPE :(h + 0.5*magSqr(U)) - p
);




Info << "\nReading Arbitrary Mesh Object\n";
pointVectorField& MDN_ = MDN;
volScalarField& rho_theo_ = rho_theo;
volVectorField& U_theo_ = U_theo;
volScalarField& p_theo_ = p_theo;
volScalarField& rho_ = rho;
volVectorField& rhoU_ = rhoU;
volScalarField& rhoE_ = rhoE;
volScalarField& JW_ = JW;
arbMesh arbitraryMesh(MDN_, mesh, runTime, rho_theo_, U_theo_, p_theo_, rho_, rhoU_, rhoE_, JW_);
arbMesh& aMsh = arbitraryMesh;




#include "updateFields.H"

forAll(C, cell)
{
  // update rhoU and rhoE
  rhoU[cell] = rho[cell]*U[cell];
  rhoE[cell] = rho[cell]*h[cell] -p[cell];
}





// Create numeric flux
Info << "\nInstanciation of flux\n";
autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New(p, U, T, thermo() );
basicNumericFlux& dbnsFlux = dbnsFluxPtr();
// Create mass flux alias for easier coupling with other code components
const surfaceScalarField& phi = dbnsFlux.rhoFlux();



  



Info << "\nWrite-up of some variables\n" << endl;
rho.write();
U.write();
p.write();
MDN.write();
rho_TALE.write();
U_TALE.write();
