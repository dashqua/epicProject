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
  const pointMesh& pMesh_;
  const fvMesh& mesh_;
  //const label& face_;
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
        explicit arbMesh(pointVectorField& MDN, const fvMesh& mesh, const Time& runTime)
	:
	  MDN_(MDN),
	  pMesh_(MDN.mesh()),
	  mesh_(mesh)
        {
	  // for deltaw and Uwn, it is necessary to keep initial areas
	  
	}

        void hello();
        void createFields(pointVectorField MDN);
        void updateMeshDisplacement(const scalar t); 
        pointVectorField MDN();
        scalar jw(vector& xyz);
        scalar deltaw(vector&  xyz, label& face);
        scalar Uwn(vector& xyz, vector& n);
        scalar Shift(scalar lambda, vector& xyzOwn, vector& xyzNei);
        scalar getMagSf(vectorList x);
        scalarList cross(scalarList A, scalarList B);
        vector vw(vector& xyz);
        //- Construct from components
        //arbMesh(const scalar& data);

        //- Construct from Istream
        //arbMesh(Istream&);

        //- Construct as copy
        //arbMesh(const arbMesh&);

	  /*
	  pointVectorField xN_
	  (
	     IOobject("xN_", mesh_),
	     pMesh_,
	     dimensionedVector("xN_", dimensionSet(0,1,0,0,0,0,0), vector::zero)
	  );
	  */

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
