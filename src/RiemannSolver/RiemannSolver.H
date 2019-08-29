/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RiemannSolver

Description

SourceFiles
    RiemannSolverI.H
    RiemannSolver.C
    RiemannSolverIO.C

\*---------------------------------------------------------------------------*/

#ifndef RiemannSolver_H
#define RiemannSolver_H

//#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// Forward declaration of friend functions and operators
class RiemannSolver;


/*---------------------------------------------------------------------------*\
                         Class RiemannSolver Declaration
\*---------------------------------------------------------------------------*/

class RiemannSolver
:
    public baseClassName
{
    // Private data

        //- Description of data_
  //dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
  //RiemannSolver(const RiemannSolver&);

        //- Disallow default bitwise assignment
  //void operator=(const RiemannSolver&);


public:

    // Static data members

        //- Static data staticData
  //static const dataType staticData;


    // Constructors

        //- Construct null
        RiemannSolver();

        //- Construct from components
  //RiemannSolver(const dataType& data);

        //- Construct from Istream
  //RiemannSolver(Istream&);

        //- Construct as copy
  //RiemannSolver(const RiemannSolver&);


    // Selectors

        //- Select null constructed
        static autoPtr<RiemannSolver> New();


    //- Destructor
    ~RiemannSolver();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

  //void operator=(const RiemannSolver&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

  //friend Istream& operator>>(Istream&, RiemannSolver&);
  //friend Ostream& operator<<(Ostream&, const RiemannSolver&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "RiemannSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
