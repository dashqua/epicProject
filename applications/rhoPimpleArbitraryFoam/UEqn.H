//Update arbitrary velocity
forAll(mesh.cells(), cell)
{
  scalar x = mesh.C[cell].x();
  scalar y = mesh.C[cell].y();
  Uarb[cell] = vector
    (
     4.*pi/T *Foam::sin(pi*x/10.)*Foam::sin(pi*y*2./15.)*Foam::cos(2.*pi*t/T),
     6*.pi/T *Foam::sin(pi*x/10.)*Foam::sin(pi*y*2./15.)*Foam::cos(4.*pi*t/T),
     0
    );
}

// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
 fvm::ddt(rho, U) + fvm::div(phi, U) - fvm::div(phi, Uarb)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
