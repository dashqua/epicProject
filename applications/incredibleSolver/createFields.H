Info << "Reading thermophysical properties\n" << endl;
autoPtr<basicPsiThermo> thermo
(
 basicPsiThermo::New(mesh)
);

Info << "\n Reading mesh information\n";
Info << " [+] Reading cell centers\n";
const surfaceVectorField& Sf = mesh.Sf();
//const volVectorField& C = mesh.C();
Info << " [+] Reading pMesh\n";
pointMesh pMesh(mesh);



// Primitive variables
Info << "\nReading primitive variables\n";
Info << " [+] Reading temperature T\n"; const volScalarField& T = thermo->T();
Info << " [+] Reading pressure p\n"   ; volScalarField& p = thermo->p(); 
Info << " [+] Reading density rho\n"  ;
volScalarField rho
(
 IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 thermo->rho()
);
Info << " [+] Reading velocity U\n"   ;
volVectorField U
(
 IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
 mesh
 );
Info << " [+] Reading enthalpy h and internal energy e\n";
volScalarField& he = thermo->h(); // he();
volScalarField& h  = he;
volScalarField& small_e  = he;
if (he.name() == "e")
  {
    volScalarField& small_e = he;
    volScalarField h_obj
      (
       IOobject("h", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
       small_e + p/rho
       );
    volScalarField& h = h_obj; Info << h.name();  //  just for remove some warnings
  }
 else // (he.name() == "h")
   {
     volScalarField& h = he;
     volScalarField e_obj
       (
	IOobject("e", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
	h - p/rho
	);
     volScalarField& small_e = e_obj; Info << small_e.name();  //  just for remove some warnings
   }
Info << " [+] Reading Total Energy E\n";
volScalarField E
(
 IOobject("E", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 small_e + 0.5*magSqr(U)
);
Info << " [+] Reading Cv\n";
volScalarField Cv = thermo->Cv();


// check dim U, E and p
///Info << U.dimensions;
//U.write(); E.write(); p.write();


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
 rho * (h + 0.5*magSqr(U)) - p //rho*E
);





Info << "\nReading arbitrary Mesh variables\n";
Info << " [+] Reading Mesh Displacement Vector\n";
pointVectorField MDN
(
 IOobject("MDN", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 pMesh,
 dimensionedVector("MDN", dimLength, vector::zero)
);

// write mesh velocity field on mesh not pmesh for now
Info << " [+] Reading arbitrary velocity Utilde\n";
volVectorField U_tilde
(
 IOobject("U_tilde", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
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







// Creation of homemade numerical flux
Info << "\nReading numerical flux\n";
Info << " [+] Reading rho Flux\n";
surfaceScalarField rhoFlux
(
 IOobject("rhoFlux", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 linearInterpolate( rho*U ) & Sf
);
Info << " [+] Reading rhoU Flux\n";
surfaceVectorField rhoUFlux
(
 IOobject("rhoUFlux", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 rhoFlux * linearInterpolate( U )
);
Info << " [+] Reading rhoE Flux\n";
surfaceScalarField rhoEFlux
(
 IOobject("rhoEFlux", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
 rhoFlux * linearInterpolate(Cv*T + magSqr(U))//linearInterpolate( U * (E+p) ) & Sf
);






// Create Gradient fields
Info << "\nReading Gradient Fields\n";
Info << " [+] Reading grad P\n"; volVectorField gradP(fvc::grad(p));
Info << " [+] Reading grad U\n"; volTensorField gradU(fvc::grad(U));
Info << " [+] Reading grad T\n"; volVectorField gradT(fvc::grad(T));


// Create numeric flux
Info << "\nInstanciation of DBNS flux\n";
autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New(p, U, T, thermo() );
basicNumericFlux& dbnsFlux = dbnsFluxPtr();



Info << "\nReading Arbitrary Mesh Object\n";
pointVectorField& MDN_    = MDN;
volScalarField& rho_      = rho;
volVectorField& U_        = U;
volScalarField& p_        = p;
volScalarField& E_        = E;
volScalarField& h_        = h;
volScalarField& rho_theo_ = rho_theo;
volVectorField& U_theo_   = U_theo;
volScalarField& p_theo_   = p_theo;
volVectorField& rhoU_     = rhoU;
volScalarField& rhoE_     = rhoE;
volVectorField& U_tilde_  = U_tilde;
surfaceScalarField& rhoFlux_  = rhoFlux;//dbnsFlux.rhoFlux();
surfaceVectorField& rhoUFlux_ = rhoUFlux;//dbnsFlux.rhoUFlux();
surfaceScalarField& rhoEFlux_ = rhoEFlux;//dbnsFlux.rhoEFlux();
arbMesh arbitraryMesh(MDN_, thermo(), rho_, U_, E_, p_, h_, rho_theo_, U_theo_, p_theo_, rhoU_, rhoE_, U_tilde_, rhoFlux_, rhoUFlux_, rhoEFlux_);
arbMesh& aMsh = arbitraryMesh;


// Create mass flux alias for easier coupling with other code components
const surfaceScalarField& phi = rhoFlux;

// Now that every variables is read, they are corrected to comply with the initial conditions of arbMesh
//aMsh.correctInitialVariables();

// Create Riemann Solver
Info << "\nInitializing the Riemann Solver\n";
RiemannSolver RS(thermo(), U_, E_, p_, h_, T);

Info << "\nWrite-up of some variables\n" << endl;
rho.write();
U.write();
p.write();
MDN.write();
p_theo.write();
U_theo.write();
T.write();
