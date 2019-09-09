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

pointVectorField Foam::arbMesh::MDN()
{
  return MDN_;
}

scalar Foam::arbMesh::jw(vector& xyz)
{
    // Get Jw
  scalar pii = constant::mathematical::pi;
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
  scalar pii = constant::mathematical::pi;
  scalar Tper = 2;
  res[0] = ptX + 2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(2*pii*t/Tper);
  res[1] = ptY + 3/2 * Foam::sin(pii*ptX/10) * Foam::sin(2*pii*ptY/15) * Foam::sin(4*pii*t/Tper);
  res[2] = 0;
  return res;
}


scalar Foam::arbMesh::Shift(scalar lambda, vector& xyzOwn, vector& xyzNei, label& face, vector Sf)
{
  // MAJ: APPLY EVERY METHODS TO AVERAGE OF XYZOWN AND XYZNEI
  vector xyzAve = (xyzOwn + xyzNei) / 2;
  scalar magSf_ = mag(Sf);
  vector n = Sf/magSf_;
  scalar deltaw = this->deltaw(magSf_, face);
  scalar jw = this->jw(xyzAve);
  scalar Uwn = (this->vw(xyzOwn) + this->vw(xyzNei))/2 & n ;
  scalar lambda_res = (deltaw/jw) * (lambda - Uwn);
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

vector arbMesh::phix(vector& xyz)
{
  scalar pii = constant::mathematical::pi;
  scalar T = 2;
  scalar t = mesh_.time().value();
  scalar x = xyz[0];
  scalar y = xyz[1];
  return vector( \
		x + 2*Foam::sin(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/T),   \
		y + 3/2*Foam::sin(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(4*pii*t/T), \
		0);
}

vector arbMesh::vw(vector& xyz)
{
  scalar pii = constant::mathematical::pi;
  scalar T = 2;
  scalar t = mesh_.time().value();
  scalar x = xyz[0];
  scalar y = xyz[1];
  return vector( \
    4*pii/T * Foam::sin(pii*x/10) * Foam::sin(2*pii*y/15) * Foam::cos(2*pii*t/T), \
    6*pii/T * Foam::sin(pii*x/10) * Foam::sin(2*pii*y/15) * Foam::cos(4*pii*t/T), \
    0);
}

vector arbMesh::c(label& cell)
{
  
  scalar x = this->mesh_.C()[cell].x();
  scalar y = this->mesh_.C()[cell].y();
  scalar z = this->mesh_.C()[cell].z();
  vector pos = vector(x,y,z);
  vector& posref = pos;
  vector c =this->U_theo_[cell] - this->vw(posref);
  return c;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void arbMesh::updateFields()
{
  Info << "\nUpdating Fields\n" << endl;
  Info << " [+] Updating TALE variables\n";
  Info << " [+] Updating U, rho and p fields\n";
  scalar M  = 0.5;
  scalar M2 = M*M;
  scalar I  = 5.0;
  scalar I2 = I*I;
  scalar r  = 1.5;
  scalar theta = Foam::atan(0.5);
  const volScalarField Cv = this->Cv_;
  const volScalarField Cp = this->Cp_;
  scalar Rhoinf = 1;
  scalar Uinf = 0.8944;
  scalar Vinf = 0.4472;
  scalar Pinf = 3;
  scalar pii  = constant::mathematical::pi;
  scalar pii2 = pii*pii;
  scalar t = mesh_.time().value();
  volVectorField C = mesh_.C();
  forAll(C, cell)
    {
      scalar gamma = Cp[cell] / Cv[cell];
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      //scalar x3 = C[cell].z();
      //vector pos = vector(x1, x2, x3);
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);

      // Updating theoretical U, rho and p    
      rho_theo_[cell] = Rhoinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/(r*r)) , 1/(gamma-1));
      U_theo_[cell] = vector(
	     Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/(r*r))/2 ),
	     Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/(r*r))/2 ),
	     0);
      p_theo_[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/(r*r)) , gamma/(gamma-1));
    }
  Info << " [+] Update Mesh Displacement\n";
  forAll(mesh_.points(), ptI)
    {
      scalar pX  = mesh_.points()[ptI].x();
      scalar pY  = mesh_.points()[ptI].y();
      scalar pZ  = mesh_.points()[ptI].z();
      vector pos = vector(pX, pY, pZ);
      // Updating Mesh Displacement Vector
      MDN_[ptI] = this->phix(pos);    
    }  
}

tensor arbMesh::Fw(vector& xyz)
{
  scalar t = mesh_.time().value();
  scalar x = xyz[0], y = xyz[1];//, z = xyz[2];
  scalar pii = constant::mathematical::pi, Tper = 2;
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

scalar arbMesh::detFw(vector& xyz)
{
  tensor Fw = this->Fw(xyz);
  return( Fw.xx() * Fw.yy() - Fw.xy() * Fw.yx());
}

tensor arbMesh::transposeFw(vector& xyz)
{
  tensor Fw = this->Fw(xyz);
  tensor transposeFw(Fw.xx(), Fw.yx(), 0,
		     Fw.xy(), Fw.yy(), 0,
		     0,0,1);
  return transposeFw;
}

tensor arbMesh::invFw(vector& xyz)
{
  tensor Fw = this->Fw(xyz);
  tensor invFw(  Fw.yy() , -Fw.xy(), 0,
		 -Fw.yx(), Fw.xx() , 0,
		       0 ,       0 , 1	   );
  return invFw / this->detFw(xyz);
}

tensor arbMesh::Hw(vector& xyz)
{
  tensor Fw = this->Fw(xyz);
  tensor Hw(  Fw.yy() , -Fw.yx(), 0,
	      -Fw.xy(), Fw.xx() , 0,
	      0 ,       0 , 1	   );	    
  return Hw;
}

tensor arbMesh::transposeHw(vector& xyz)
{
  tensor Fw = this->Fw(xyz);
  tensor Hw(  Fw.yy() , -Fw.xy(), 0,
	      -Fw.yx(), Fw.xx() , 0,
	      0 ,       0 , 1	   );	    
  return Hw;
}

tensor arbMesh::inverseTransposeHw(vector&xyz)
{
  return 1/this->detFw(xyz) * Fw(xyz);
}

/////////////////////////////////////////////////////////////////////////////

void arbMesh::computeTALEfromEUL()
{
  // If not precised, variables are assumed to be EULERIAN
  // e.g. rhoU is rho_EUL * U_EUL
  volVectorField C = mesh_.C();

  /*
  forAll(mesh_.owner(), face)
    {
      scalar x = mesh_.Cf()[face].x();
      scalar y = mesh_.Cf()[face].y();
      scalar z = mesh_.Cf()[face].z();
      vector pos = vector(x,y,z);

      // Update Fluxes      /!\ similar pattern here
      rhoFlux_[face]  = this->detFw(pos) * this->transposeHw(pos) & rhoFlux_[face];
      rhoUFlux_[face] = this->detFw(pos) * this->transposeHw(pos) & rhoUFlux_[face];
      rhoEFlux_[face] = this->detFw(pos) * this->transposeHw(pos) & rhoEFlux_[face];
    }



  /*
  volVectorField tempRhoFlux = rho_ * U_;
  volVectorField tempRhoUFlux = U_; //its OK its temp
  volVectorField tempRhoEFlux = U_; // Cv_ * T_ + 0.5*magSqr(U_);
  forAll(C, cell)
    {
      scalar x = C[cell].x();
      scalar y = C[cell].y();
      scalar z = C[cell].z();
      vector pos = vector(x,y,z);

      // Update Primitive Fields
      rho_[cell] = this->detFw(pos) * rho_[cell];
      U_[cell] = this->detFw(pos) * U_[cell];
      E_[cell] = this->detFw(pos) * E_[cell];
      
      // Update Conserved Fields
      rhoU_[cell] = rho_[cell]*U_[cell];
      rhoE_[cell] = rho_[cell]*E_[cell];

      // Update Fluxes via temporary variables
      tempRhoFlux[cell]  = this->detFw(pos) * this->transposeHw(pos) & ( rho_[cell] * U_[cell] );
    }
  rhoFlux_  = (linearInterpolate(tempRhoFlux) & mesh_.Sf());
  rhoUFlux_ = rhoFlux_ * linearInterpolate(U_);
  rhoEFlux_ = rhoFlux_ * linearInterpolate(Cv_*T_ + 0.5*magSqr(U_));

  */
}


void arbMesh::computeEULfromTALE()
{

}

void arbMesh::correctInitialVariables()
{
  Info << "\nCorrection of Initial Variables\n" << endl;
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
  scalar t = mesh_.time().value();
  volVectorField C = mesh_.C();
  Info << " [+] Correcting theoretical and conserved variables\n";
  forAll(C, cell)
    {
      scalar gamma = Cp_[cell]/Cv_[cell];
      scalar x1 = C[cell].x();
      scalar x2 = C[cell].y();
      //scalar x3 = C[cell].z();
      //vector pos = vector(x1, x2, x3);
      scalar v1 = Uinf*Foam::cos(theta);
      scalar v2 = Vinf*Foam::sin(theta);
      Info << Cp_[cell] << "   " << Cv_[cell];//pow( 1 - I2*M2*(gamma-1)/(8*pii2)  , 1./(gamma-1) ) ; Info << "top";      
      // Correcting theoretical U, rho and p    
      rho_theo_[cell] = Rhoinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , 1/(gamma-1)); Info << "1";
      U_theo_[cell] = vector(
	     Uinf * ( Foam::cos(theta)- I*(x2-v2*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
	     Vinf * ( Foam::sin(theta)- I*(x1-v1*t)/(2*pii*r) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2)/2 ),
	     0);
      p_theo_[cell] = Pinf * Foam::pow( 1 - I2*M2*(gamma-1)/(8*pii2) * Foam::exp((1-(x1-v1*t)*(x1-v1*t)-(x2-v2*t)*(x2-v2*t))/r2) , gamma/(gamma-1));

      // Correcting Primitive Fields       
      rho_[cell] = rho_theo_[cell];
      U_[cell]   = U_theo_[cell];
      p_[cell]   = p_theo_[cell]; 

      // Correcting Conserved Fields
      rhoU_[cell] = rho_[cell] * U_[cell];
      rhoE_[cell] = rho_[cell] * E_[cell]; // (h_[cell] + 0.5*magSqr(U_));

      // Fluxes don't need to be initialized because they are computed afterwards
    }
  //  /!\ AT THE MOMENT, H AND E ARE NOT CORRECTED
  //E_   = p_/(gamma-1) + 0.5 * rho_ * magSqr(U_);
  //h_   = p_/(rho_*(gamma-1)) + 0.5 * magSqr(U_);

  Info << " [+] test: Correcting Boundary Conditions\n";
  U_.correctBoundaryConditions();
  h_.correctBoundaryConditions();
  p_.correctBoundaryConditions();
  
  Info << " [+] Correcting Mesh Displacement\n";
  forAll(mesh_.points(), ptI)
    {
      scalar pX  = mesh_.points()[ptI].x();
      scalar pY  = mesh_.points()[ptI].y();
      scalar pZ  = mesh_.points()[ptI].z();
      vector pos = vector(pX, pY, pZ);
      // Updating Mesh Displacement Vector
      MDN_[ptI] = this->phix(pos);    
    }    
}

void arbMesh::computeFluxALEfromEUL()
{
  const unallocLabelList& owner     = mesh_.owner();
  const unallocLabelList& neighbour = mesh_.neighbour();
  
    
  forAll(owner, face)
  {
    // Step 0: Decoding quantities
    const label own = owner[face];
    const label nei = neighbour[face];
    const scalar pLeft     = p_[own] ;//+ .... ;
    const vector ULeft     = U_[own] ;//+ ..... ;
    const scalar TLeft     = T_[own] ;//+ ... ;
    const scalar hLeft     = h_[own];
    const scalar RLeft     = Cp_[own] - Cv_[own];
    const scalar CvLeft    = Cv_[own];
    const scalar CpLeft    = Cp_[own];
    const scalar gammaLeft = CpLeft/CvLeft;
    //
    const scalar pRight     = p_[nei] ;//+ .... ;
    const vector URight     = U_[nei] ;//+ ..... ;
    const scalar TRight     = T_[nei] ;//+ ... ;
    const scalar hRight     = h_[nei];
    const scalar RRight     = Cp_[nei] - Cv_[nei];
    const scalar CvRight    = Cv_[nei];
    const scalar CpRight    = Cp_[nei];
    const scalar gammaRight = CpLeft/CvLeft;    
    //     code right values here before compiling /!\

    const scalar rhoLeft     = pLeft/(RLeft*TLeft);
    const scalar rhoRight    = pRight/(RRight*TRight);
    const scalar rhoLeftSqr  = Foam::sqrt(max(rhoLeft,SMALL));
    const scalar rhoRightSqr = Foam::sqrt(max(rhoRight,SMALL));
    const scalar eLeft       = CvLeft*TLeft+0.5*magSqr(ULeft);
    const scalar eRight      = CvRight*TRight+0.5*magSqr(URight);
    const scalar wLeft       = rhoLeft/(rhoLeftSqr+rhoRightSqr);
    const scalar wRight      = rhoRight/(rhoLeftSqr+rhoRightSqr);
    const scalar kappaLeft   = Cp_[own]/Cv_[own];
    const scalar kappaRight  = Cp_[nei]/Cv_[nei];
      //
    const scalar rhoTilde    = Foam::sqrt(max( rhoLeft*rhoRight ,SMALL));
    const vector UTilde      = ULeft*wLeft + URight*wRight;
    const scalar hTilde      = hLeft*wLeft + hRight*wRight;
    const scalar gammaTilde  = gammaLeft*wLeft + gammaRight*wRight;
    const scalar qTildeSquare= mag(UTilde);
    const scalar kappaTilde  = kappaLeft*wLeft + kappaRight*wRight;
    
    const vector Sf = mesh_.Sf()[face];
    const scalar magSf = mesh_.magSf()[face];
    const vector normalVector = Sf / magSf;
    const scalar contrVTilde = UTilde & normalVector;
    const scalar contrVLeft  = ULeft & normalVector;
    const scalar contrVRight = URight & normalVector;

    const scalar deltaP = pRight - pLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = deltaU & normalVector;
    const scalar deltaRho = rhoRight - rhoLeft;
    const scalar aTilde = sqrt(max( (kappaTilde-1)*(hTilde-0.5*qTildeSquare) ,SMALL));
    
    // Step 1: Computation of eigenvalues
    vector centerLeft_ = mesh_.Cf()[own]; vector& centerLeft = centerLeft_;
    vector centerRight_= mesh_.Cf()[nei]; vector& centerRight = centerRight_;
    const vector arbULeft   = this->vw(centerLeft);    
    const vector arbURight  = this->vw(centerRight);
    const vector arbUTilde  = arbULeft*wLeft + arbURight*wRight;  
    const scalar lambda_1   = 2*Foam::sqrt(max( rhoTilde,SMALL)) * contrVTilde;
    const scalar lambda_234 = lambda_1 / 2;
    const scalar lambda_5   = lambda_234 / gammaTilde;
    //const scalar lambdaMax  = max(max(lambda_1,lambda_234),lambda_5);

    // Step 2: Computation of eigenvectors
    const scalar K1_1   = 1;
    const vector K1_234 = UTilde;
    const scalar K1_5   = (hTilde+ (gammaTilde-1)*qTildeSquare) / (2*gammaTilde -1);
    //
    const scalar K23_1   = 0;
    const vector K23_234 = vector(1,1,1) - normalVector;
    const scalar K23_5   = UTilde & ( vector(1,1,1) - normalVector);
    //
    const scalar K4_1    = 0;
    const vector K4_234  = UTilde;
    const scalar K4_5    = qTildeSquare;
    //
    const scalar K5_1    = 0;
    const vector K5_234  = vector(0,0,0);
    const scalar K5_5    = 1;

    // Step 2b: some alpha *eigenv
    const scalar alphaK23_1   = 0;
    const vector alphaK23_234 = rhoTilde * ( deltaU - deltaContrV*normalVector  );
    const scalar alphaK23_5   = rhoTilde * (( UTilde & deltaU ) - contrVTilde*deltaContrV);
    
    // Step 3: Recalculation of wave strengths
    // ATM : they are the same as for eulerian
    const scalar alpha1   = (deltaP - rhoTilde*aTilde*deltaContrV)/(2.0*sqr(aTilde));
    //const scalar alpha23= rhoTilde * ((deltaU & vector(1,1,1))  - deltaContrV);
    const scalar alpha4   = deltaRho - deltaP/sqr(aTilde);
    const scalar alpha5   = (deltaP + rhoTilde*aTilde*deltaContrV)/(2.0*sqr(aTilde));

    // Step 4: Assembly of flux differences
    const scalar DF1_1    = lambda_1 * alpha1 * K1_1;
    const vector DF1_234  = lambda_1 * alpha1 * K1_234;
    const scalar DF1_5    = lambda_1 * alpha1 * K1_5;
    //
    const scalar DF23_1   = lambda_234 * alphaK23_1;
    const vector DF23_234 = lambda_234 * alphaK23_234;
    const scalar DF23_5   = lambda_234 * alphaK23_5;
    //
    const scalar DF4_1    = lambda_234 * alpha4 * K4_1;
    const vector DF4_234  = lambda_234 * alpha4 * K4_234;
    const scalar DF4_5    = lambda_234 * alpha4 * K4_5;
    //
    const scalar DF5_1    = lambda_5 * alpha5 * K5_1;
    const vector DF5_234  = lambda_5 * alpha5 * K5_234;
    const scalar DF5_5    = lambda_5 * alpha5 * K5_5;

    // Step 5: Update of fluxes
    const scalar fluxALE1   = 0.5 * (DF1_1 + DF23_1 + DF4_1 + DF5_1);
    const vector fluxALE234 = 0.5 * (DF1_234 + DF23_234 + DF4_234 + DF5_234);
    const scalar fluxALE5   = 0.5 * (DF1_5 + DF23_5 + DF4_5 + DF5_5);
    
    // Step 5b: Apply the updates on flux
    rhoFlux_  -= fluxALE1 * magSf;
    rhoUFlux_ -= fluxALE234 * magSf;
    rhoEFlux_ -= fluxALE5 * magSf;
  }
}


void arbMesh::updateCoordinates()
{
  
}

/*
surfaceScalarField arbMesh::FluxTALEfromEUL(const surfaceScalarField& rhoFlux)
{
  // Getting TALE from EUL Scalar flux <=> multiply by JH^-T = F^-1
  surfaceScalarField rhoFlux_TALE = rhoFlux;
  forAll(mesh_.Cf(), face)
    {
      vector pos = mesh_.Cf()[face];
      rhoFlux_TALE[face] = 
    }
}
  */

//surfaceVectorField arbMesh::FluxTALEfromEUL(const surfaceVectorField& rhoUFlux)
//{


//}




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
