pointMesh pMesh(mesh);
const pointField meshpoints = mesh.points();

pointVectorField tabulated_dat
(
 IOobject("tabulated_dat",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
 pMesh,
 dimensionedVector("tabulated_dat", dimLength, vector::zero)
);
tabulated_dat.write();

volVectorField U
(
 IOobject("U",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
 mesh,
 dimensionedVector("U", dimLength, vector::zero)
 );
U.write();
