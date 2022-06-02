/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::LESModels::sstTransitionDES

Group
    grpDESTurbulence

Description
    k-omega-SST DES turbulence model for incompressible and compressible flows

    Reference:
    \verbatim
        Strelets, M. (2001)
        Detached Eddy Simulation of Massively Separated Flows,
        39th AIAA Aerospace Sciences Meeting and Exhibit, Reno, NV
    \endverbatim

Note
    The default values of the DES constants implemented are code-specific
    values calibrated for OpenFOAM using decaying isotropic turbulence, and
    hence deviate slightly from the values suggested in the reference
    publication.

SourceFiles
    sstTransitionDES.C

\*---------------------------------------------------------------------------*/

#ifndef sstTransitionDES_H
#define sstTransitionDES_H

#include "DESModel.H"
#include "sstTransitionBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        class sstTransitionDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class sstTransitionDES
:
    public sstTransitionBase<DESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        sstTransitionDES(const sstTransitionDES&) = delete;

        //- No copy assignment
        void operator=(const sstTransitionDES&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar CDESkom_;
            dimensionedScalar CDESkeps_;


    // Protected Member Functions

        //- Blending for CDES parameter
        virtual tmp<volScalarField> CDES(const volScalarField& F1) const
        {
            return this->blend(F1, CDESkom_, CDESkeps_);
        }

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;

        //- Return epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("sstTransitionDES");


    // Constructors

        //- Construct from components
        sstTransitionDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~sstTransitionDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "sstTransitionDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
