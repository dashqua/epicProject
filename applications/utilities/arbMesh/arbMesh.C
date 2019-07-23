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

void Foam::arbMesh::updateMeshDisplacement(const scalar t)
{
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
      dA = mag(this->getMagSf(coords));
      coords[0] = this->apply_mapping(coords[0]);
      coords[1] = this->apply_mapping(coords[1]);
      coords[2] = this->apply_mapping(coords[2]);
      coords[3] = this->apply_mapping(coords[3]);
      da = mag(this->getMagSf(coords));
      //Info << "size4:\n da=" << da << "\ndA=" << dA << "\nmagSf_=" << magSf_<<endl;
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
      dA = mag(this->getMagSf( coords ));
      coords[0] = this->apply_mapping(coords[0]);
      coords[1] = this->apply_mapping(coords[1]);
      coords[2] = this->apply_mapping(coords[2]);
      da = mag(this->getMagSf( coords ));
      //Info << "mesh.face:     " << mesh_.faces()[face] << "\nface: " << face << endl;
      //Info << "size3:\n da=" << da << "\n dA=" << dA
      //<< "\n mesh_.magSf()=" << mesh_.magSf()[face] << "\n magSf_=" << magSf_
      //  <<endl<<endl<<endl;
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

scalar Foam::arbMesh::Uwn(vector& xyzOwn, vector& xyzNei, label& face) //previous cchi
{
  // /!\ IMPLEMENTED ONLY FOR TRIANGLE ATM
  const label& node0 = mesh_.faces()[face][0];
  const label& node1 = mesh_.faces()[face][1];
  const label& node2 = mesh_.faces()[face][2];
  vectorList coords(3);
  coords[0] = mesh_.points()[node0];
  coords[1] = mesh_.points()[node1];
  coords[2] = mesh_.points()[node2];
  vector Sf = this->getMagSf( coords );
  scalar magSf = mag(Sf);
  vector n = Sf/magSf;

  vector xyzAve = (xyzOwn + xyzNei) / 2 ;
  vector vw = this->vw(xyzAve);
  return Foam::dot( vw , n );
}

scalar Foam::arbMesh::Shift(scalar lambda, vector& xyzOwn, vector& xyzNei, label& face, vector Sf)
{
  // MAJ: APPLY EVERY METHODS TO AVERAGE OF XYZOWN AND XYZNEI
  vector xyzAve = (xyzOwn + xyzNei) / 2;
  scalar magSf_ = mag(Sf);
  //vector n = Sf/magSf_;
  scalar deltaw = this->deltaw(magSf_, face), \
    jw = this->jw(xyzAve),                    \
    Uwn = this->Uwn(xyzOwn, xyzNei, face),    \
    lambda_res = (deltaw/jw) * (lambda - Uwn);
  /*
  //left
  scalar lambdaLeft = (deltaw/jw) * (lambda - Uwn); 
  //right
  jw = this->jw(xyzNei);
  Uwn = this->Uwn(xyzNei, n);
  scalar lambdaRight = (deltaw/jw) * (lambda - Uwn);
  //stab
  scalar lambdaStab = 0;
  */
  return lambda_res;
}

vector Foam::arbMesh::getMagSf(vectorList x)
{ 
  vectorList iso(x.size());
  vector tanVecZeta = vector::zero;
  vector tanVecEta = vector::zero;
  //scalar a = 0;
  vector n = vector::zero;
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
      n = (tanVecZeta ^ tanVecEta)/mag(tanVecZeta ^ tanVecEta); //Sf
      //a = 4.0*mag(tanVecZeta ^ tanVecEta); //magSf
    }
  if (x.size()==3)
    {
      /*
      vector A1 = x[1]-x[0], A2 = x[2]-x[0]; 
      a = Foam::sqrt( (A2.y()*A1.z()-A2.y()*A1.z()) * (A2.y()*A1.z()-A2.y()*A1.z()) +\
		      (A1.x()*A2.z()-A2.x()*A1.z()) * (A1.x()*A2.z()-A2.x()*A1.z()) +\
		      (A1.x()*A2.y()-A2.x()*A1.y()) * (A1.x()*A2.y()-A2.x()*A1.y()) \
		      ) / 2 ;
      */
      vector BA = x[1] - x[0], CA = x[2] - x[0];
      n = (BA ^ CA) / 2 ;
      //a = mag(n);
    }
  return n;
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void arbMesh::updateFields(scalar t)
{
  this->updateMeshDisplacement(t);
  //const scalar Tper = 2;
  const scalar pii = Foam::mathematicalConstant::pi;
  // mean velocity field;
  vector v_moy = vector::zero; v_moy[0] = 0.8944; v_moy[1] = 0.4472;
    // analytics fields update;
  const scalar I = 5.0, r = 1.5, gamma = 1.330, M = sqrt(v_moy[0]*v_moy[0]+v_moy[1]*v_moy[1]),
    theta = 26.56*pii/180, x10 = 0, x20 = 0, v1 = .8944, v2 = .4472;
  forAll(this->mesh_.C(),cell)
  {
    scalar x1 = this->mesh_.C()[cell].x();
    scalar x2 = this->mesh_.C()[cell].y();

    vector tmp_vec = vector(x1, x2, 0);
    JW_[cell] = this->jw( tmp_vec );
    rho_theo_[cell] = 1 * Foam::pow( 1 - I*I*(gamma-1)*M*M/(8*pii*pii) * Foam::exp(f_fun(x1,x2,t,v_moy)) , 1/(gamma-1) );
    p_theo_[cell] = 3 * Foam::pow( 1 - I*I*(gamma-1)*M*M/(8*pii*pii) * Foam::exp(f_fun(x1,x2,t,v_moy)) , gamma/(gamma-1));
    // INCORRECT U THEO
    U_theo_[cell] = vector(
			   sqrt(v_moy[0]*v_moy[0]+v_moy[1]*v_moy[1]) * (Foam::cos(theta)- I*(x2-x20-v2*t)/(2*pii*r)*Foam::exp(f_fun(x1,x2,t,v_moy)/2)),
			   sqrt(v_moy[0]*v_moy[0]+v_moy[1]*v_moy[1]) * (Foam::sin(theta)+ I*(x1-x10-v1*t)/(2*pii*r)*Foam::exp(f_fun(x1,x2,t,v_moy)/2)),
			   0 );
  }  
}

scalar arbMesh::f_fun(scalar x1,scalar x2, scalar t, vector v)
{
  scalar r = 1.5, x10 = 0, x20 = 0, v1 = v[0], v2 = v[1];
  return (1 - (x1-x10-v1*t)*(x1-x10-v1*t) - (x2-x20-v2*t)*(x2-x20-v2*t) )/(r*r);
}

tensor arbMesh::Fw(scalar x, scalar y, scalar z, scalar t)
{
  scalar pii = Foam::mathematicalConstant::pi, Tper = 2;
  tensor F(
	   1+pii/5*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/Tper),
	   4*pii/15*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(2*pii*t/Tper),
	   0,
	   3*pii/20*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(4*pii*t/Tper),
	   1+pii/5*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(4*pii*t/Tper),
	   0,
	   0,
	   0,
	   1
  );
  return F;
}

tensor arbMesh::cross(tensor& A, tensor &B)
{
  tensor C = A;
  for (int i=0; i<A.dim1(); i++)
    {
      for (int j=0; j<B.dim2(); j++)
	{
	  C[i,j]
	}
    }
}

/*
tensor arbMesh::Hw(scalar x, scalar y, scalar z, scalar t)
{
  return this->Fw(x,y,z,t)
}*/

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
