// Update coordinates
// Info << "test: " << xN << endl;

//Info << typeof(MDN) << endl << endl << endl << endl << endl << endl << endl;

const scalar pii = Foam::mathematicalConstant::pi;
const scalar Tper = 2;
const scalar t = runTime.value();
const pointField& meshpoints = mesh.points();
forAll(meshpoints, ptI)
{
  scalar ptX = meshpoints[ptI].x();
  scalar ptY = meshpoints[ptI].y();
  MDN[ptI] = vector(							\
     ptX + 2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(2*pii*t/Tper),   \
     ptY + 3/2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(4*pii*t/Tper), \
     0                                                                                      \
  );
}

// Info << meshpoints[0].x() << endl;

uN = xN - XN;
//uN.write();

if (runTime.outputTime())
  {
    MDN.write();
  }
