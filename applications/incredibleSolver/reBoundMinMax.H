void Foam::reBoundMinMax
(
 volScalarField& vsf,
 dimensionedScalar& vsf0,
 dimensionedScalar& vsf1
)
{
  scalar minVsf = min(vsf).value();
  scalar maxVsf = max(vsf).value();

  if (minVsf < vsf0.value() || maxVsf > vsf1.value())
  {
    Info << "bounding " << vsf.name()
	 << ", min: " << minVsf
	 << " max: " << maxVsf
	 << " average: " << gAverage(vsf.internalField())
	 << " average prim: " << gAverage(vsf.primitiveField())
	 << endl;
  }

  if (minVsf < vsf0.value())
  {
    vsf.primitiveFieldRef() = max
      (
       max
       (
	vsf.primitiveFieldRef(),
	fvc::average( max(vsf,vsf0) )().primitiveField() *
	pos( vsf0.value()- vsf.primitiveField() )
       ),
        vsf0.value()
      );

    Info << "new min: " << gMin(vsf.internalField()) << endl;
    vsf.correctBoundaryConditions();
    vsf.boundaryFieldRef() = max(vsf.boundaryField(), vsf0.value());
  }

  if (maxVsf > vsf1.value())
  {
    vsf.primitiveFieldRef() = min
    (
      min
      (
       vsf.primitiveFieldRef(),
       fvc::average(min(vsf,vsf1))().primitiveField()
       *neg(vsf1.value() - vsf.primitiveField())
       //
       + pos(vsf1.value() - vsf.primitiveField())*vsf1.value()
      )
    );
    Info << "new max: " << gMax(vsf.internalField()) << endl;
    vsf.correctBoundaryConditions();
    vsf.boundaryFieldRef() = min(vsf.boundaryField(), vsf1.value());
  }
}
