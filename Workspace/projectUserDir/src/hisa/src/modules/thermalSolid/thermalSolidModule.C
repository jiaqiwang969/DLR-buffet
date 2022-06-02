/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2021 Oliver Oxtoby

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

\*---------------------------------------------------------------------------*/

#include "thermalSolidModule.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(thermalSolidModule, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void thermalSolidModule::updateAniAlpha(const tmp<volVectorField>& tkappaByCp)
{
    aniAlpha_->primitiveFieldRef() =
#if OPENFOAM >= 1812
        coordinates_->transformPrincipal(mesh_->C(), tkappaByCp());
#else
        coordinates_->R().transformVector(tkappaByCp());
#endif
    forAll(aniAlpha_->boundaryField(), patchi)
    {
        aniAlpha_->boundaryFieldRef()[patchi] =
#if OPENFOAM >= 1812
            coordinates_->transformPrincipal
            (
                mesh_->C().boundaryField()[patchi],
                tkappaByCp->boundaryField()[patchi]
            );
#else
            aniAlpha_->boundaryField()[patchi].patchInternalField();
#endif
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalSolidModule::thermalSolidModule
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    solverModule(name),
    time_(t),
    dict_(dict)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalSolidModule::initialise()
{
    const Time& runTime = time_;
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    // Detect steady-state analysis
    const dictionary& ddtControls = mesh.schemesDict().subDict("ddtSchemes");
    wordList ddtToc (ddtControls.toc());
    steadyState_ = false;
    forAll(ddtToc,s)
    {
        word ddtScheme(ddtToc[s]);
        word ddtSchemeLastWord;
        const tokenList& tokens = ddtControls.lookup(ddtToc[s]);
        if (tokens.last().isWord() && tokens.last().wordToken() == "steadyState")
        {
            if (ddtToc[s] == "default" || ddtToc[s] == "ddt(betav*rho,h)")
            {
                steadyState_ = true;
            }
        }
    }
    if (steadyState_)
    {
        Info << "Steady-state analysis detected" << nl << endl;
        solnControl_.set(new simpleControl(mesh));
    }
    else
    {
        Info << "Transient analysis detected" << nl << endl;
        solnControl_.set(new pimpleControl(mesh));
    }
}

scalar thermalSolidModule::timeStepScaling(const scalar& maxCoNum)
{
    #include "solidRegionDiffusionNo.H"

    return maxDi_/DiNum;
}

void thermalSolidModule::beginTimeStep()
{
    maxDi_ = time_.controlDict().lookupOrDefault<scalar>("maxDi", GREAT);
    moveMeshOuterCorrectors_ = 
        solnControl_().dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", true);
    firstOuterCorrector_ = true;
}

void thermalSolidModule::outerIteration()
{
    #include "setSolidFields.H"

    if (isA<dynamicFvMesh>(mesh) && (firstOuterCorrector_ || moveMeshOuterCorrectors_))
    {
        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

        dynamicCast<dynamicFvMesh&>(mesh).update();

        if (mesh.changing())
        {
            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;
        }
    }

    #include "solveSolid.H"

    firstOuterCorrector_ = false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    addToRunTimeSelectionTable
    (
        solverModule,
        thermalSolidModule,
        dictionary
    );
}

// ************************************************************************* //

