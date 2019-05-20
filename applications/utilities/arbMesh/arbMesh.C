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
  Info << "mesh.face.size: " << mesh_.faces()[0].size() <<endl;
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

scalar Foam::arbMesh::jw(vector& xyz)
{
    // Get Jw
    scalar pii = Foam::mathematicalConstant::pi;
    scalar Tper = 2;
    scalar t = mesh_.time().value();
    scalar x = xyz[0];
    scalar y = xyz[1];
    //scalar z = xyz[2];    
    scalar jw = ( 1 + pii/5*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/Tper)  )    * \
      ( 1 + pii/5*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(4*pii*t/Tper)  )              - \
      ( 3*pii/20*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(4*pii*t/Tper)   )              * \
      ( 4*pii/15*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(2*pii*t/Tper)   )              ;
    return jw;
}

scalar Foam::arbMesh::deltaw(scalar magSf_, label& face)
{
  //scalar id = mesh_.findCell(point(x,y,z));
  scalar dA = 0, da = 0;//mesh_.magSf()[id];
  //forAll(mesh_.faces(), face)
  //{
  if (mesh_.faces()[face].size()==4)
    {
      const label& node0 = mesh_.faces()[face][0];
      const label& node1 = mesh_.faces()[face][1];
      const label& node2 = mesh_.faces()[face][2];
      const label& node3 = mesh_.faces()[face][3];
      vectorList coords(4);
      coords[0] = mesh_.points()[node0]; // coord0;
      coords[1] = mesh_.points()[node1]; // coord1;
      coords[2] = mesh_.points()[node2]; //coord2;
      coords[3] = mesh_.points()[node3]; //coord3;
      dA = this->getMagSf(coords);
      coords[0] = this->apply_mapping(coords[0]);
      coords[1] = this->apply_mapping(coords[1]);
      coords[2] = this->apply_mapping(coords[2]);
      coords[3] = this->apply_mapping(coords[3]);
      da = this->getMagSf(coords);
      Info << "size4:\n da=" << da << "\ndA=" << dA << "\nmagSf_=" << magSf_<<endl;
      return da/dA;      
    }
  else if (mesh_.faces()[face].size()==3)
    {
      const label& node0 = mesh_.faces()[face][0];
      const label& node1 = mesh_.faces()[face][1];
      const label& node2 = mesh_.faces()[face][2];
      vectorList coords(3);
      coords[0] = mesh_.points()[node0]; //coord0; //this->apply_mapping(coord0);
      coords[1] = mesh_.points()[node1]; //coord1; //this->apply_mapping(coord1);
      coords[2] = mesh_.points()[node2]; // mesh_.points()[node0]; //coord2; //this->apply_mapping(coord2);
      dA = this->getMagSf( coords );
      coords[0] = this->apply_mapping(coords[0]);
      coords[1] = this->apply_mapping(coords[1]);
      coords[2] = this->apply_mapping(coords[2]);
      da = this->getMagSf( coords );
      Info << "mesh.face:     " << mesh_.faces()[face] << endl;
      Info << "size3:\n da=" << da << "\n dA=" << dA
	   << "\n mesh_.magSf()=" << mesh_.magSf()[face] << "\n magSf_=" << magSf_
	   <<endl<<endl<<endl;
      /*
      vectorList t_test(4);
      t_test[0] = vector(0,0,0);
      t_test[1] = vector(1,0,0);
      t_test[2] = vector(1,1,0);
      t_test[3] = vector(0,1,0);
      scalar test = this->getMagSf( t_test );
      Info << "~~~~~  " << test << endl <<endl;
      */
      return da/dA;
    }
  else
    {
      Info << "Problem in deltaw; area is null" << endl;
      return 1; // 0;
    }
}

vector Foam::arbMesh::apply_mapping(vector coord)
{
  scalar t = mesh_.time().value();
  vector res = coord;
  scalar ptX = res[0];  // apply MDN_ displacement for the point coord
  scalar ptY = res[1];
  scalar pii = Foam::mathematicalConstant::pi;
  scalar Tper = 2;
  res[0] = ptX + 2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(2*pii*t/Tper);
  res[1] = ptY + 3/2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(4*pii*t/Tper);
  res[2] = 0;
  return res;
}

scalar Foam::arbMesh::Uwn(vector& xyz, vector n) //previous cchi
{
  vector vw = this->vw(xyz);
  return Foam::dot( vw , n );
}

scalar Foam::arbMesh::Shift(scalar lambda, vector& xyzOwn, vector& xyzNei, label& face, vector Sf)
{
  //    return (lambdaw/jw) * (lambda - cchi);
  // following remark is initially for flux::::
  // the physical flux is sum of the left and right numerical fluxes
  // in addition to a correction term (stabilization)
  // /!\ DELTAW IS THE SAME FOR THE TWO SIDES
  scalar deltaw = 0, jw = 0, Uwn = 0;
  scalar magSf_ = mag(Sf);
  deltaw = this->deltaw(magSf_, face);
  vector n = Sf/magSf_;
  //left
  jw = this->jw(xyzOwn);
  Uwn = this->Uwn(xyzOwn, n);
  scalar lambdaLeft = (deltaw/jw) * (lambda - Uwn); 
  //right
  jw = this->jw(xyzNei);
  Uwn = this->Uwn(xyzNei, n);
  scalar lambdaRight = (deltaw/jw) * (lambda - Uwn);
  //stab
  scalar lambdaStab = 0;
  return (lambdaLeft + lambdaRight)/2 + lambdaStab;
}

scalar Foam::arbMesh::getMagSf(vectorList x)
{ 
  vectorList iso(x.size());
  vector tanVecZeta = vector::zero;
  vector tanVecEta = vector::zero;
  scalar a = 0;
  if (x.size()==4)
    {
      iso[0] = vector(-1,-1,0);
      iso[1] = vector(1,-1,0);
      iso[2] = vector(1,1,0);
      iso[3] = vector(-1,1,0);
    
      for (int i =0; i<4; i++)
	{
	  tanVecZeta += x[i] * (iso[i].x() / 4.0);
	  tanVecEta  += x[i] * (iso[i].y() / 4.0);
	}
      //vector n = (tanVecZeta ^ tanVecEta)/mag(tanVecZeta ^ tanVecEta); //Sf
      a = 4.0*mag(tanVecZeta ^ tanVecEta); //magSf
    }
  if (x.size()==3)
    {
      /*
      vectorList vlist(3);
      vlist[0] = x[0] - x[1];//x[0];
      vlist[1] = x[1] - x[2];//x[1];
      vlist[2] = x[2] - x[0];//x[2];
      vectorList& vlist_ref = vlist;
      Info << "before : " << vlist_ref << endl;
      sortVlist(vlist_ref);
      Info << "after : " << vlist_ref << endl;
      scalar cc = mag(vlist_ref[2]), bb = mag(vlist_ref[1]), aa = mag(vlist_ref[0]);
      Info << aa << endl << bb << endl << cc << endl;
      Info << "square value : " << (aa+(bb+cc)) * (cc-(aa-bb)) * (cc+(aa-bb)) * (aa+(bb-cc)) << endl;
      return Foam::sqrt( mag( (aa+(bb+cc)) * (cc-(aa-bb)) * (cc+(aa-bb)) * (aa+(bb-cc)) ) );
      */
      /*
      Info << "C11:  " << ( (x[1].x()-x[0].x())*(x[2].y()-x[0].y()) - (x[2].x()-x[0].x())*(x[1].y()-x[0].y()) ) / 2
	   << "  C22:" << ( (x[2].y()-x[0].y())*(x[1].x()-x[0].x()) - (x[1].y()-x[0].y())*(x[2].x()-x[0].x()) ) / 2
	   << endl;
      */
      a =  ( (x[1].x()-x[0].x())*(x[2].y()-x[0].y()) - (x[2].x()-x[0].x())*(x[1].y()-x[0].y()) ) / 2;
    }
  return a;
}

void Foam::arbMesh::sortVlist(vectorList& vlist)
{
  if (mag(vlist[0]) > mag(vlist[1]))
    {
      vector& tmp = vlist[0];
      vlist[0] = vlist[1];
      vlist[1] = tmp;
    }
  if ( mag(vlist[0]) > mag(vlist[2]) ) 
    {
      vector& tmp = vlist[0];
      vlist[0] = vlist[2];
      vlist[2] = tmp;
    }
  if (    (mag(vlist[2]) < mag(vlist[1])) && (mag(vlist[2]) > mag(vlist[0]))  )
    {
      vector& tmp = vlist[1];
      vlist[1] = vlist[2];
      vlist[2] = tmp;
    }
}

scalar arbMesh::facto(scalar n)
{
  if (n>1)
    {
      return n*this->facto(n-1);
    }
  else
    {
      return 1;
    }
}

scalar arbMesh::arccos(scalar x)
{
  scalar pii = Foam::mathematicalConstant::pi;
  /*
  scalar a = pii/2;
  for (int i=0; i<n; i++)
    {
      a -= this->facto(2*n)/pow( this->facto(n)*Foam::pow(2,n) , 2 ) * pow(x, 2*n +1)/(2*n +1);
      } */
  return pii/2 - (x + x*x*x/6 + 3*x*x*x*x*x/40 + 15*x*x*x*x*x*x*x/336);
}

vector Foam::arbMesh::vec(vector A, vector B)
{
  return A - B;
}

scalar Foam::arbMesh::len_vec(vector A, vector B)
{  // returns length between point A and B
  vector AB = this->vec(A, B);
  scalar x = AB.x(), y = AB.y();
  return Foam::sqrt( x*x + y*y );
}

scalarList Foam::arbMesh::cross(scalarList A, scalarList B)
{
  scalarList res(3);
  res[0] = (A[1]*B[2] - A[2]*B[1]);
  res[1] = (A[2]*B[0] - A[0]*B[2]);
  res[2] = (A[0]*B[1] - A[1]*B[0]);
  return res;
}

vector Foam::arbMesh::vw(vector& xyz)
{
  scalar pii = Foam::mathematicalConstant::pi;
  scalar Tper = 2;
  scalar t = mesh_.time().value();
  scalar x = xyz[0];
  scalar y = xyz[1];
  return vector( \
    4*pii/Tper * Foam::sin(pii*x/10) * Foam::sin(2*pii*y/15) * Foam::cos(2*pii*t/Tper), \
    6*pii/Tper * Foam::sin(pii*x/10) * Foam::sin(2*pii*y/15) * Foam::cos(4*pii*t/Tper), \
    0);
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
