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

Class
    Foam::arbMesh

Description

SourceFiles
    arbMeshI.H
    arbMesh.C
    arbMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef arbMesh_H
#define arbMesh_H

#include "polyMesh.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "fvCFD.H"
#include "pointMesh.H"
#include "GeometricField.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "vectorList.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class volMesh;
  
/*---------------------------------------------------------------------------*\
                          Class arbMesh Declaration
\*---------------------------------------------------------------------------*/

class arbMesh
//:
//    public polyMesh
{
private:
  pointVectorField& MDN_;
  basicThermo& thermo_;
  volScalarField& rho_;
  volVectorField& U_;
  volScalarField& p_;
  volScalarField& E_;
  const volScalarField& T_;
  volScalarField& h_;
  volScalarField& rho_theo_;
  volVectorField& U_theo_;
  volScalarField& p_theo_;
  volVectorField& rhoU_;
  volScalarField& rhoE_;
  volVectorField& U_tilde_;
  const fvMesh& mesh_;
  const volScalarField& Cv_;
  const volScalarField& Cp_;
  surfaceScalarField& rhoFlux_;
  surfaceVectorField& rhoUFlux_;
  surfaceScalarField& rhoEFlux_;


  
    // Private data
        //- Description of data_
        //scalar data_;
    // Private Member Functions
        //- Disallow default bitwise copy construct
        //arbMesh(const arbMesh&);
        //- Disallow default bitwise assignment
        //void operator=(const arbMesh&);
public:

  // Static data members
  //- Static data staticData
  //static const dataType staticData;
  //- Construct from Field
  explicit arbMesh(pointVectorField& MDN,
		   basicThermo& thermo,
		   volScalarField& rho,		   
		   volVectorField& U,
		   volScalarField& E,
		   volScalarField& p,
		   volScalarField& h,
		   volScalarField& rho_theo,
		   volVectorField& U_theo,
		   volScalarField& p_theo,
		   volVectorField& rhoU,
		   volScalarField& rhoE,
		   volVectorField& U_tilde,
		   surfaceScalarField& rhoFlux,
		   surfaceVectorField& rhoUFlux,
		   surfaceScalarField& rhoEFlux )
    :
    MDN_(MDN),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    p_(p),
    E_(E),
    T_(thermo.T()),
    h_(h),
    rho_theo_(rho_theo),
    U_theo_(U_theo),
    p_theo_(p_theo),
    rhoU_(rhoU),
    rhoE_(rhoE),
    U_tilde_(U_tilde),
    rhoFlux_(rhoFlux),
    rhoUFlux_(rhoUFlux),
    rhoEFlux_(rhoEFlux),
    mesh_(U.mesh()),
    Cv_(thermo_.Cv()),
    Cp_(thermo_.Cp())
    { }
  //rhoFlux_  = rhoFlux;
    //rhoUFlux_ = rhoUFlux;
    //rhoEFlux_ = rhoEFlux;
    //mesh_ = p.mesh();
  
  ////////////////////////////////////////////////
  
  pointVectorField MDN();
  scalar jw(vector& xyz);
  scalar deltaw(scalar msgSf_, label& face);
  vector apply_mapping(vector coord);
  scalar Shift(scalar lambda, vector& xyzOwn, vector& xyzNei, label& face, vector n);
  vector getMagSf(vectorList x);


  ////////////////////////////////////////////////
  vector phix(vector& xyz);                     // Arbitrary Mapping
  vector vw(vector& xyz);                       // Velocity of the mapping
  vector c(label& cell);                        // convective velocity
  
  ////////////////////////////////////////////////
  void updateFields();                          // Update theoretical, TALE and MDN fields
  tensor Fw(vector& xyz);                       // Deformation Gradient
  scalar detFw(vector& xyz);                    // Jacobian of Def Grad (Vol Map)
  tensor transposeFw(vector& xyz);              // transpose of Def Grad 
  tensor invFw(vector& xyz);                    // inverse of Def Grad
  tensor Hw(vector& xyz);                       // Area Map JF^-T
  tensor transposeHw(vector& xyz);              // transpose of Area Map
  tensor inverseTransposeHw(vector& xyz);       // inverseTranspose of Area Map
  
  ////////////////////////////////////////////////
  void computeTALEfromEUL();                    // Compute TALE variables from EUL variables
  void computeEULfromTALE();                    // Compute EUL variables from TALE variables
  void correctInitialVariables();               // Correct Variables to comply with initial cond


  //////////////////////////////////////////////// Compute ALE Fluxes - substract Delta tilde F_i
  void computeFluxALEfromEUL();

  
  ////////////////////////////////////////////////
  //surfaceScalarField FluxTALEfromEUL(const surfaceScalarField& rhoFlux);
  //surfaceVectorField FluxTALEfromEUL(const surfaceVectorField& rhoUFlux);
  
  //- Construct from components
        //arbMesh(const scalar& data);

        //- Construct from Istream
        //arbMesh(Istream&);

        //- Construct as copy
        //arbMesh(const arbMesh&);

    // Selectors

        //- Select null constructed
        static autoPtr<arbMesh> New();


    //- Destructor
    //~arbMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

  //void operator=(const arbMesh&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

  //friend Istream& operator>>(Istream&, arbMesh&);
  //friend Ostream& operator<<(Ostream&, const arbMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "arbMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
