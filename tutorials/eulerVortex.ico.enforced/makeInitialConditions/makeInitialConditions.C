#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //                                                                                                                                                                                                                                                                                                              

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "initContinuityErrs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //                                                                                                                                                                                                                                                                                                              


  Info << " [+] Correcting theoretical variables\n";
  Info << " [+] Correcting U, rho and p fields\n";
  scalar M  = 0.5;
  scalar M2 = M*M;
  scalar I  = 5.0;
  scalar I2 = I*I;
  scalar r  = 1.5;
  scalar r2 = r*r;
  scalar theta = Foam::atan(0.5);
  scalar Rhoinf = 1;
  scalar Uinf = 0.8944;
  scalar Vinf = 0.4472;
  scalar Pinf = 3.;
  scalar pii  = constant::mathematical::pi;
  scalar pii2 = pii*pii;
  scalar t = mesh.time().value();
  volVectorField C = mesh.C();


  Info << " [+] Instanciation of boundary conditions for rho, p, U, T\n";
  volScalarField rho
  (
   IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("rho", dimensionSet(1,-3,0,0,0,0,0), 1)
   );

  volScalarField p
  (
   IOobject("p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("p", dimensionSet(0,2,-2,0,0,0,0), 3)
   );

  volVectorField U
  (
   IOobject("U", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE), mesh,
   dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),Foam::vector(.8944,.4472,0))
  );

  volScalarField T
  (
   IOobject("T", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh,
   dimensionedScalar("T", dimensionSet(0,0,0,1,0,0,0), 300)
  );


  forAll(C, cell)
    {
      scalar gamma = 1.3;//Cp_[cell]/Cv_[cell];                                                                                                                                                                                               
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);

      // Correcting theoretical U, rho and p                                                                                                                                                                                                  
      rho[cell] = Rhoinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , 1/(gamma-1));
      U[cell] = vector(
             Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
             0);
      p[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , gamma/(gamma-1)) ;
    }

  
  U.write();
  p.write();
  rho.write();
  T.write();


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //                                                                                                                                                                                                                                                                                                              


    Info<< "End\n" << endl;

    return 0;
}
