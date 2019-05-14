/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "arbMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

/*
Foam::arbMesh Foam::arbMesh(const fvMesh& mesh)
{
  Info << "hello arbMesh" << endl;
}
*/



void Foam::arbMesh::hello()
{
  Info << "hello" << endl;
}

void Foam::arbMesh::createFields(pointVectorField MDN_)
{
  Info << endl;
}

void Foam::arbMesh::updateMeshDisplacement(const scalar t)
{
  Info << endl;
  const scalar Tper = 2;
  const scalar pii = Foam::mathematicalConstant::pi;
  const pointField& meshpoints = mesh_.points();
  forAll(meshpoints,ptI)
  {
    scalar ptX = meshpoints[ptI].x();
    scalar ptY = meshpoints[ptI].y();
    MDN_[ptI] = vector(							\
       ptX + 2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(2*pii*t/Tper),   \
       ptY + 3/2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(4*pii*t/Tper), \
       0                                                                                      \
    );

  }
}

pointVectorField Foam::arbMesh::MDN()
{
  return MDN_;
}

scalar Foam::arbMesh::jw(scalarList& xyz)
{
    // Get Jw
    scalar pii = Foam::mathematicalConstant::pi;
    scalar Tper = 2;
    scalar t = mesh_.time().value();
    scalar x = xyz[0];
    scalar y = xyz[1];
    scalar z = xyz[2];    
    scalar jw = ( 1 + pii/5*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/Tper)  ) * \
      ( 1 + pii/5*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(4*pii*t/Tper)  )              - \
      ( 3*pii/20*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(4*pii*t/Tper)   )              * \
      ( 4*pii/15*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(2*pii*t/Tper)   )              ;
    return jw;
}

scalar Foam::arbMesh::deltaw(scalarList& xyz)
{
  scalar x = xyz[0];
  scalar y = xyz[1];
  scalar z = xyz[2];
  scalar id = mesh_.findCell(point(x,y,z));
  scalar dA = mesh_.magSf()[id];
  
}

scalar Foam::arbMesh::Uwn(scalarList& xyz) //previous cchi
{

}

scalar Foam::arbMesh::Shift(scalar lambda, scalarList& xyzOwn, scalarList& xyzNei)
{
  //    return (lambdaw/jw) * (lambda - cchi);
  // TEMPORARY
  scalar x = xyzOwn[0];
  scalar y = xyzOwn[1];    
  //return ( this->deltaw(x, y)/this->jw(x, y) ) * (lambda - this->Uwn(x, y)) ;
}

scalar Foam::arbMesh::getMagSf(scalarList& A, scalarList& B, scalarList& C)
{
  scalarList Sf = this->getSf(A, B, C);
  return Foam::sqrt(Sf[0]*Sf[0] + Sf[1]*Sf[1] + Sf[2]*Sf[2]);
}

scalarList Foam::arbMesh::getSf(scalarList& A, scalarList& B, scalarList& C)
{
  scalarList BA(3), CA(3);
  BA[0] = B[0] - A[0];
  BA[1] = B[1] - A[1];
  BA[2] = B[2] - A[2];
  CA[0] = C[0] - A[0];
  CA[1] = C[1] - A[1];
  CA[2] = C[2] - A[2];
  return this->cross(BA, CA);
}

scalarList Foam::arbMesh::cross(scalarList A, scalarList B)
{
  scalarList res(3);
  res[0] = (A[1]*B[2] - A[2]*B[1]);
  res[1] = (A[2]*B[0] - A[0]*B[2]);
  res[2] = (A[0]*B[1] - A[1]*B[0]);
  return res;
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
Foam::arbMesh::arbMesh()
{
  data_ = 0;
}


Foam::arbMesh::arbMesh(const scalar& data)
{
  data_ = data;
}


Foam::arbMesh::arbMesh(const arbMesh& am)
{
  data_ = am.data();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::arbMesh> Foam::arbMesh::New()
{
    return autoPtr<arbMesh>(new arbMesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::arbMesh::~arbMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

scalar Foam::arbMesh::data()
{
  return data_;
}

void Foam::arbMesh::operator=(const arbMesh& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::arbMesh::operator=(const Foam::arbMesh&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

*/

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
