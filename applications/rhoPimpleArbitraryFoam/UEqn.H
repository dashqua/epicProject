// Solve the Momentum equation
MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
 fvm::ddt(rho, U) + fvm::div(phi, U)//fvm::div(phi, convVel) //fvm::div(phi, U) - fvm::div(phi, Uarb)
  + MRF.DDt(rho, U)
 //+ turbulence->divDevRhoReff(U)
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
