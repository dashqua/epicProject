pointMesh pMesh(mesh);
const pointField meshpoints = mesh.points();

pointVectorField tabulated_dat
(
 IOobject("tabulated_dat",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
 pMesh,
 dimensionedVector("tabulated_dat", dimLength, vector::zero)
);
tabulated_dat.write();

volVectorField U
(
 IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
 mesh
 );

