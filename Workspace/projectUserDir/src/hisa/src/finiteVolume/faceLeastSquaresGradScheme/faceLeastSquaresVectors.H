/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
    Copyright (C) 2011 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::faceLeastSquaresVectors

Description
    Least-squares gradient scheme vectors

SourceFiles
    faceLeastSquaresVectors.C

Authors
    Oliver Oxtoby
    Johan Heyns
        Council for Scientific and Industrial Research, South Africa

    Based on leastSquares library of OpenFOAM

\*---------------------------------------------------------------------------*/

#ifndef faceLeastSquaresVectors_H
#define faceLeastSquaresVectors_H

#ifdef BLUECFD
#include "MeshObject.T.H"
#else
#include "MeshObject.H"
#endif
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faceLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class faceLeastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, faceLeastSquaresVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("faceLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit faceLeastSquaresVectors(const fvMesh&);


    //- Destructor
    virtual ~faceLeastSquaresVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;


        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
