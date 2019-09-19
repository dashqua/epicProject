/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dynamicArbitraryFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicArbitraryFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicArbitraryFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicArbitraryFvMesh::dynamicArbitraryFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).optionalSubDict(typeName + "Coeffs")
    ),
    amplitude_(readScalar(dynamicMeshCoeffs_.lookup("amplitude"))),
    frequency_(readScalar(dynamicMeshCoeffs_.lookup("frequency"))),
    refPlaneX_(readScalar(dynamicMeshCoeffs_.lookup("refPlaneX"))),
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    Info<< "Performing a dynamic mesh calculation: " << endl
        << "amplitude: " << amplitude_
        << " frequency: " << frequency_
        << " refPlaneX: " << refPlaneX_ << endl;

    initialPoints_ = stationaryPoints_;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicArbitraryFvMesh::~dynamicArbitraryFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicArbitraryFvMesh::update()
{
  /*
    scalar scalingFunction =
        0.5*
        (
            ::cos(constant::mathematical::twoPi*frequency_*time().value())
          - 1.0
        );
  */
    Info<< "Mesh scaling. Time = " << time().value() << " scaling: "
      //<< scalingFunction << endl;
	<< "arbitrary scaling" << endl;


    scalar pi = constant::mathematical::pi;
    scalar T = 2.;
    scalar t = this->time().value();
    
    pointField newPoints = stationaryPoints_;    
    forAll(newPoints, point)
      {
	//Info << newPoints[point] << endl ;
	scalar x = initialPoints_[point].x();
	scalar y = initialPoints_[point].y();
	scalar z = initialPoints_[point].z();
	newPoints[point] = vector
	  (
	   x + 2*Foam::sin(pi*x/10.)*Foam::sin(2*pi*y/15.)*Foam::sin(2*pi*t/T),
	   y + 3./2. *Foam::sin(pi*x/10.)*Foam::sin(2*pi*y/15.)*Foam::sin(4*pi*t/T),
	   0
	   );
      }
    
    fvMesh::movePoints(newPoints);

    lookupObjectRef<volVectorField>("U").correctBoundaryConditions();

    return true;
}


// ************************************************************************* //
