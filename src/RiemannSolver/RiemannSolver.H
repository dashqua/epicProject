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

#include "polyMesh.H"
#include "GeometricField.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "volPointInterpolation.H"
#include "basicThermo.H"
#include "MDLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RiemannSolver Declaration
\*---------------------------------------------------------------------------*/

class RiemannSolver
{
private:

  const fvMesh& mesh_;
  basicThermo& thermo_;
  volVectorField& U_;
  volScalarField& E_;
  volScalarField& p_;
  volScalarField& h_;
  const volScalarField& T_;
  
public:
  // Constructor
  explicit RiemannSolver
  (
    basicThermo& thermo,
    volVectorField& U,
    volScalarField& E,
    volScalarField& p,
    volScalarField& h,
    const volScalarField& T
   )
    :
    mesh_(U.mesh()),
    thermo_(thermo),
    U_(U),
    E_(E),
    p_(p),
    h_(h),
    T_(T)
  {}

  // Destructor
  ~RiemannSolver()
  {}

  // Methods
  void computeFlux(surfaceScalarField& , surfaceVectorField& , surfaceScalarField&, volVectorField&, volTensorField&, volVectorField&);

  void evaluateFluxInternal(scalar&, vector&, scalar&, const scalar&, const scalar&, const vector&, const vector&, const scalar&, const scalar&, const scalar&, const scalar&, const scalar&, const scalar&, const vector&, const scalar&);
  
  void evaluateFluxBoundary(scalar&, vector&, scalar&, const scalar&, const scalar&, const vector&, const vector&, const scalar&, const scalar&, const scalar&, const scalar&, const scalar&, const scalar&, const vector&, const scalar&);
  
  //- Select null constructed
  //static autoPtr<RiemannSolver> New();


  













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
