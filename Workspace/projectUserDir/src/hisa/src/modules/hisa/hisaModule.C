/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa

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

#include "hisaModule.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "preconditioner.H"
#include "solver.H"
#include "fvcSmooth.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "dynamicRefineFvMesh.H"
#include "processorFvPatchField.H"
#include "wallDist.H"
#include "jacobian.H"
#include "characteristicWallPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(hisaModule, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


void hisaModule::createPreconditioners
(
    PtrList<preconditioner<2,1>>& preconditioners,
    PtrList<jacobian>& jacobians,
    const dictionary& parentDict
)
{
    if (parentDict.found("preconditioner"))
    {
        word preconditionerName(parentDict.lookup("preconditioner"));
        label nj = jacobians.size();
        label np = preconditioners.size();

        if (parentDict.found(preconditionerName))
        {
            const dictionary& dict = parentDict.subDict(preconditionerName);

            // Create preconditioner-specific Jacobian if specified
            if (dict.found("inviscidJacobian"))
            {
                jacobians.append
                (
                    new jacobian
                    (
                        parentDict.subDict(preconditionerName),
                        mesh(),
                        scalarVars_[0],
                        vectorVars_[0],
                        scalarVars_[1],
                        pThermo_(),
                        U_(),
                        ddtCoeff_(),
                        inviscid_,
#if FOUNDATION >= 8
                        turbulence_,
                        thermophysicalTransport_
#else
                        turbulence_
#endif
                    )
                );
                nj = jacobians.size();
            }

            // Recursively create further preconditioners if applicable
            createPreconditioners(preconditioners, jacobians, dict);
        }

        // Create this preconditioner; pass newly created Jacobian if created,
        // else previous one, and pointer to child preconditioner if one
        // was created otherwise NULL.
        preconditioners.append
        (
            preconditioner<2,1>::New
            (
                parentDict,
                jacobians[nj-1].matrix(),
                np > 0 ? preconditioners(np-1) : NULL
            )
        );
    }
}

template<class FieldType, class Type>
void hisaModule::parallelSyncFields(const wordList& fields)
{
    forAll(fields, i)
    {
        FieldType& f = const_cast<FieldType&>(mesh_->lookupObject<FieldType>(fields[i]));
        forAll(f.boundaryField(), patchI)
        {
            if (isA<processorFvPatchField<Type>>(f.boundaryField()[patchI]))
            {
                f.boundaryFieldRef()[patchI].initEvaluate();
            }
        }
        forAll(f.boundaryField(), patchI)
        {
            if (isA<processorFvPatchField<Type>>(f.boundaryField()[patchI]))
            {
                f.boundaryFieldRef()[patchI].evaluate();
            }
        }
    }
}

void hisaModule::redistributePar()
{
    fvMesh& mesh = mesh_();

    // Check load balance
    scalarList procLoad(Pstream::nProcs(), 0.0);
    procLoad[Pstream::myProcNo()] = mesh.nCells();

    reduce(procLoad, sumOp<List<scalar>>());

    scalar overallLoad = sum(procLoad);
    scalar averageLoad = overallLoad/double(Pstream::nProcs());

    bool balanced = true;
    for (int i = 0; i < Pstream::nProcs(); i++)
    {
        if (Foam::mag(procLoad[i] - averageLoad)/averageLoad > maxLoadImbalance_)
        {
            balanced = false;
        }
    }
    Info << "Max load imbalance: " << max(Foam::mag(procLoad-averageLoad)/averageLoad)*100.0 << " %" << endl;

    // Redistribute
    if (!balanced)
    {
        if (rebalance_)
        {
            Info << "Redistributing parallel decomposition" << endl;
            scalar timeBeforeDist = mesh.time().elapsedCpuTime();

            labelList decomposition(mesh.nCells(), 0);

            {
                IOdictionary decomposeParDict
                (
                    IOobject
                    (
                        "decomposeParDict",
                        mesh.time().system(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                // For convenience, try parallel scotch
                if (word(decomposeParDict.lookup("method")) == "scotch")
                {
                    decomposeParDict.set("method", "ptscotch");
                }
                autoPtr<decompositionMethod> decomposer
                (
                    decompositionMethod::New
                    (
                        decomposeParDict
                    )
                );
                decomposition = decomposer().decompose(mesh, mesh.cellCentres());
            }

            // Demand-driven data is cleared by fvMeshDistribute, but turbulence model has
            // stored a reference to the wall distance field. We preserve the object
            // (prevent a dangling reference) by checking it out of the object registry
            // and checking back in afterwards. The wall distance field itself stays
            // checked in and is redistributed together with all the other fields.
            wallDist* pwd = 0;
            if (mesh.foundObject<wallDist>(wallDist::typeName))
            {
                pwd = &const_cast<wallDist&>(mesh.lookupObject<wallDist>(wallDist::typeName));
                pwd->release();
            }

            {
                Info << "Distributing mesh..." << endl;

                // Create mesh re-distribution engine
                #if OPENFOAM >= 2106 or FOUNDATION >= 9
                fvMeshDistribute distributor(mesh);
                #else
                const scalar tolDim = 1e-6*mesh.bounds().mag();
                fvMeshDistribute distributor(mesh, tolDim);
                #endif

                // Re-distribute the mesh
                autoPtr<mapDistributePolyMesh> map = distributor.distribute(decomposition);

                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    // Update the refinement
                    dynamicRefineFvMesh& refineMesh = dynamic_cast<dynamicRefineFvMesh&>(mesh);
                    const_cast<hexRef8&>(refineMesh.meshCutter()).distribute(map);
                    // Update protected-cell list
#if OPENFOAM >= 1712
                    bitSet& pc = refineMesh.protectedCell();
#else
                    PackedBoolList& pc = refineMesh.protectedCell();
#endif
                    List<bool> upc(pc.size());
                    forAll(pc, cellI)
                    {
                        upc[cellI] = pc[cellI];
                    }
                    map->distributeCellData(upc);
#if OPENFOAM >= 1712
                    pc = bitSet(upc);
#else
                    pc = PackedBoolList(upc);
#endif
                }
            }

            if (pwd)
            {
                pwd->store();
            }

            // The distributor zeros the processor patches, so re-sync
            const wordList vsfs(mesh.names(volScalarField::typeName));
            const wordList vvfs(mesh.names(volVectorField::typeName));
            const wordList vtfs(mesh.names(volTensorField::typeName));
            const wordList vstfs(mesh.names(volSphericalTensorField::typeName));
            const wordList vsphtfs(mesh.names(volSymmTensorField::typeName));
            parallelSyncFields<volScalarField, scalar>(vsfs);
            parallelSyncFields<volVectorField, vector>(vvfs);
            parallelSyncFields<volTensorField, tensor>(vtfs);
            parallelSyncFields<volSymmTensorField, symmTensor>(vstfs);
            parallelSyncFields<volSphericalTensorField, sphericalTensor>(vsphtfs);

            scalarList procLoadNew (Pstream::nProcs(), 0.0);
            procLoadNew[Pstream::myProcNo()] = mesh.nCells();

            reduce(procLoadNew, sumOp<List<scalar> >());

            scalar overallLoadNew = sum(procLoadNew);
            scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

            Info << "Max load imbalance: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;

            Info<< "Execution time for mesh redistribution = "
                << mesh.time().elapsedCpuTime() - timeBeforeDist << " s\n" << endl;
        }
        else
        {
            Info << "Parallel mesh redistribution not selected" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hisaModule::hisaModule
(
    const word& name,
    const Time& t
)
:
    solverModule(name),
    time_(t)
{
}

hisaModule::hisaModule
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    solverModule(name),
    time_(t)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hisaModule::initialise()
{

    const Time& runTime = time_;

    #include "createDynamicFvMesh.H"

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
            if (ddtToc[s] == "default" || ddtToc[s] == "rhoU")
            {
                steadyState_ = true;
            }
        }
    }
    if (steadyState_)
    {
        Info << "Steady-state analysis detected" << nl << endl;
    }
    else
    {
        Info << "Transient analysis detected" << nl << endl;
    }

    residualIO defaultPseudoTol(2, 1, residualOrdering, 1e-4);
    residualIO defaultPseudoTolRel(2, 1, residualOrdering, 0.0);

    solnControl_.set
    (
        new pseudotimeControl
        (
            mesh,
            steadyState_,
            2,
            1,
            defaultPseudoTol,
            defaultPseudoTolRel
        )
    );
    if (steadyState_)
    {
        solnControl_->setCorr(runTime.startTimeIndex());
    }

    localTimestepping_ =
        solnControl_->dict().lookupOrDefault<Switch>
        (
            "localTimestepping",
            true
        );
    if (localTimestepping_)
    {
        Info << "Local timestepping selected" << nl << endl;
        localTimesteppingBounding_ =
            solnControl_->dict().lookupOrDefault<Switch>
            (
                "localTimesteppingBounding",
                true
            );

        localTimesteppingLowerBound_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "localTimesteppingLowerBound",
                0.95
            );
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ > 0 ? localTimesteppingLowerBound_ : 0.0);
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ < 0.99 ? localTimesteppingLowerBound_ : 0.99);

        // localTimesteppingUpperBound_ =
        //     solnControl_->dict().lookupOrDefault<scalar>
        //     (
        //         "localTimesteppingUpperBound",
        //         1.5
        //     );
        // localTimesteppingUpperBound_ =
        //     (localTimesteppingUpperBound_ < 1.01 ? localTimesteppingUpperBound_ : 1.01);
    }
    else
    {
        Info << "Global timestepping selected" << nl << endl;
        localTimesteppingBounding_ = false;
    }

    #include "createFields.H"

    // User config check
    #include "config.H"

    // Read or initialise pseudo Co number
    // NOTE: It is not necessarily read before every outer iteration
    // (see resetPseudo)
    if (!localTimestepping_)
    {
        pseudoCoNum_.set
        (
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "pseudoCoNum",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                dimensionedScalar
                (
                    "pseudoCoNumInit",
                    dimless,
                    solnControl_().dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
                )
             )
        );

        Info << "Initial pseudo Courant No: " << pseudoCoNum_->value() << nl << endl;
    }
    else
    {
        pseudoCoField_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "pseudoCoField",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "pseudoCoNumInit",
                    dimless,
                    solnControl_().dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
                )
            )
        );
        pseudoCoField_();

        scalar totCells = mesh.globalData().nTotalCells();
        Info<< "Initial pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField_()).value()
            << " Mean: "
            << sum(pseudoCoField_()).value()/totCells
            << " Max: "
            << max(pseudoCoField_()).value() << nl << endl;
    }

    findDebugCell();
}

scalar hisaModule::timeStepScaling(const scalar& maxCoNum)
{

    const Time& runTime = time_;
    fvMesh& mesh = mesh_();
    psiThermo& thermo = pThermo_();
    volVectorField& U = U_();

    Info << "Mesh region: " << name() << endl;

    #include "compressibleCourantNo.H"

    return maxCoNum/stabilise(CoNum, SMALL);

}

void hisaModule::beginTimeStep()
{
    // Clear out residuals from previous time step
    if (!steadyState_)
    {
        initRes_.clear();
        prevRes_.clear();
    }

    pseudotimeControl& solnControl = solnControl_();
    moveMeshOuterCorrectors_ =
        solnControl.dict().lookupOrDefault<label>
        (
            "moveMeshOuterCorrectors",
            0
        );

    // Reset pseudo Co number before every time step
    if (solnControl.dict().lookupOrDefault<bool>("resetPseudo", false) )
    {
        if (!localTimestepping_)
        {
            pseudoCoNum_->value() =
                solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
        }
        else
        {
            pseudoCoField_() == solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
        }
    }
    pseudoCoNumMin_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);
    pseudoCoNumMax_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);
    pseudoCoNumMaxIncr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMaxIncreaseFactor",
        1.25
    );
    pseudoCoNumMinDecr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMinDecreaseFactor",
        0.1
    );

    rebalance_ =
        mesh_->time().controlDict().lookupOrDefault<Switch>("rebalance", false);
    maxLoadImbalance_ =
        mesh_->time().controlDict().lookupOrDefault<scalar>("maxLoadImbalance", 0.1);
}

void hisaModule::outerIteration()
{

    #include "ptrRefs.H"

    if
    (
        solnControl.firstIter()
     || steadyState_
     || (moveMeshOuterCorrectors_ && !(solnControl.corr() % moveMeshOuterCorrectors_))
    )
    {
        if (isA<dynamicFvMesh>(mesh))
        {
            // Do any mesh changes

            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            // Work around issues with mapping of oriented fields
            surfaceVectorField Un(phi*mesh.Sf()/sqr(mesh.magSf()));
            #if OPENFOAM >= 1712
            phi.setOriented(false);
            #endif

            dynamicCast<dynamicFvMesh&>(mesh).update();

            if (mesh.changing())
            {
                Info<< "Execution time for mesh.update() = "
                    << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                    << " s" << endl;
            }
            if (mesh.topoChanging())
            {
                redistributePar();
            }

            #if OPENFOAM >= 1712
            phi.setOriented(true);
            #endif
            phi = (Un & mesh.Sf());
        }
    }

    // Store value before solve for bounding. Don't use storePrevIter here
    // as it introduces too many problems for mesh refinement/redistribution
    volScalarField rhoPrevIter("rhoPrevIter", rho);
    bounded_ = false;

    phiUp_.set
    (
        new surfaceVectorField
        (
            IOobject
            (
                "phiUp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoU.dimensions()*U.dimensions()
        )
    );
    phiEp_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiEp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoE.dimensions()*U.dimensions()
        )
    );
    Up_.set
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Up",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*U.dimensions()
        )
    );

    surfaceVectorField& phiUp = phiUp_();
    surfaceScalarField& phiEp = phiEp_();
    surfaceVectorField& Up = Up_();

    // Interpolation of primitive fields on faces
    flux.calcFlux(phi, phiUp, phiEp, Up);

    // Update residuals
    #include "residualsUpdate.H"

    // Pseudo Courant number relaxation
    setPseudoCoNum();

    // Update pseudo time step
    setPseudoDeltaT();

    // Set based on latest rPseudoDeltaT
    ddtCoeff = max(max(fvm::ddt(rho)->diag(),fvm::ddt(rhoU)->diag()),fvm::ddt(rhoE)->diag())/mesh.V();

    #include "cellDebug.H"

    // Store previous iteration values. Don't use storePrevIter here as it
    // introduces too many problems for mesh refinement/redistribution
    scalarVarsPrevIter_.clear();
    vectorVarsPrevIter_.clear();
    forAll(scalarVars_, i)
    {
        scalarVarsPrevIter_.append(new volScalarField(scalarVars_[i].name() + "PrevIter", scalarVars_[i]));
    }
    forAll(vectorVars_, i)
    {
        vectorVarsPrevIter_.append(new volVectorField(vectorVars_[i].name() + "PrevIter", vectorVars_[i]));
    }

    // Solver
    {
        residualIO defaultSolverTol(2, 1, residualOrdering, 1e-12);

        // Storage of jacobians and preconditioners
        PtrList<jacobian> jacobians;
        PtrList<preconditioner<2,1>> preconditioners;

        const dictionary& dict = mesh.solutionDict().subDict("flowSolver");
        const word solverType(dict.lookup("solver"));

        // Create main Jacobian
        jacobians.append
        (
            new jacobian
            (
                dict.subOrEmptyDict(solverType),
                mesh,
                rho,
                rhoU,
                rhoE,
                thermo,
                U,
                ddtCoeff_(),
                inviscid_,
#if FOUNDATION >= 8
                turbulence_,
                thermophysicalTransport_
#else
                turbulence_
#endif
            )
        );

        // Recursively create needed preconditioners and their jacobians (if
        // applicable)
        createPreconditioners
        (
            preconditioners,
            jacobians,
            dict.subOrEmptyDict(solverType)
        );

        // Create solver
        autoPtr<solver<2,1>> sol =
            solver<2,1>::New
            (
                dict,
                jacobians[0].matrix(),
                preconditioners.size() ? preconditioners(0) : NULL,
                defaultSolverTol
            );

        // Call solver
        initRes_.clear();
        sol->solve
        (
            scalarVars_,
            vectorVars_,
            scalarResiduals_,
            vectorResiduals_,
            initRes_
        );
        solnControl_->setResidual(initRes_());
    }

    if (localTimesteppingBounding_)
    {
        volScalarField rhoMin =
            localTimesteppingLowerBound_*scalarVarsPrevIter_[0];
        // volScalarField rhoMax = 2*scalarVarsPrevIter_[0];
        volScalarField eMin =
            localTimesteppingLowerBound_*
            (scalarVarsPrevIter_[1]/scalarVarsPrevIter_[0]
            - 0.5*magSqr(vectorVarsPrevIter_[0]/scalarVarsPrevIter_[0]));
        // volScalarField eMax =
        //     100.0*(scalarVarsPrevIter_[1]/scalarVarsPrevIter_[0]
        //     - 0.5*magSqr(vectorVarsPrevIter_[0]/scalarVarsPrevIter_[0]));
        volScalarField eTemp(rhoE/rho - 0.5*magSqr(rhoU/rho));

        volScalarField factor
        (
            IOobject
            (
                "factor",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("factor", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(mesh.owner(), facei)
        {
            label own = mesh.owner()[facei];
            label nei = mesh.neighbour()[facei];

            if
            (
                (rho[own] < rhoMin[own]) //|| (rho[own] > rhoMax[own])
             || (eTemp[own] < eMin[own]) //|| (eTemp[own] > eMax[own])
             || (eTemp[own] < SMALL)
            )
            {
                factor[own] = min(0.5, factor[own]);
                factor[nei] = min(0.75, factor[nei]);
            }
            if
            (
                (rho[nei] < rhoMin[nei]) //|| (rho[nei] > rhoMax[nei])
             || (eTemp[nei] < eMin[nei]) //|| (eTemp[nei] > eMax[nei])
             || (eTemp[nei] < SMALL)
            )
            {
                factor[nei] = min(0.5, factor[nei]);
                factor[own] = min(0.75, factor[own]);
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                scalarField rhoNei =
                    rho.boundaryField()[patchi].patchNeighbourField();
                scalarField rhoMinNei =
                    rhoMin.boundaryField()[patchi].patchNeighbourField();
                // scalarField rhoMaxNei =
                //     rhoMax.boundaryField()[patchi].patchNeighbourField();
                scalarField eTempNei =
                    eTemp.boundaryField()[patchi].patchNeighbourField();
                scalarField eMinNei =
                    eMin.boundaryField()[patchi].patchNeighbourField();
                // scalarField eMaxNei =
                //     eMax.boundaryField()[patchi].patchNeighbourField();
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, bfacei)
                {
                    if
                    (
                        (rho[fc[bfacei]] < rhoMin[fc[bfacei]])
                    //  || (rho[fc[bfacei]] > rhoMax[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < eMin[fc[bfacei]])
                    //  || (eTemp[fc[bfacei]] > eMax[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.5, factor[fc[bfacei]]);
                    }
                    if
                    (
                        (rhoNei[bfacei] < rhoMinNei[bfacei])
                    //  || (rhoNei[bfacei] > rhoMaxNei[bfacei])
                     || (eTempNei[bfacei] < eMinNei[bfacei])
                    //  || (eTempNei[bfacei] > eMaxNei[bfacei])
                     || (eTempNei[bfacei] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.75, factor[fc[bfacei]]);
                    }
                }
            }
        }

        pseudoCoField_() *= factor;
    }

    phiUp_.clear();
    phiEp_.clear();
    Up_.clear();

    // Update fields
    #include "updateFields.H"

    // Correct turbulence
    turbulence_->correct();
#if FOUNDATION >= 8
    thermophysicalTransport_->correct();
#endif
}

void hisaModule::findDebugCell()
{
    // Code adapted from findRefCell.C
    const dictionary& dict = mesh_->time().controlDict();
    if (dict.found("debugCell"))
    {
        cellDebugging_ = true;

        if (!Pstream::parRun() || Pstream::myProcNo() == readLabel(dict.lookup("debugProcessor")))
        {
            debugCell_ = readLabel(dict.lookup("debugCell"));
            if (debugCell_ < 0 || debugCell_ >= mesh_->nCells())
            {
                debugCell_ = -1;
                Pout << "Debug cell " << debugCell_ << " is out of range." << endl;
            }
        }
        else
        {
            debugCell_ = -1;
        }

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        if (returnReduce<label>(hasRef, sumOp<label>()) != 1)
        {
            WarningInFunction
                << "debugCell not found." << endl;
            cellDebugging_ = false;
        }
    }
    else if (dict.found("debugPoint"))
    {
        point debugPoint(dict.lookup("debugPoint"));

        // Try fast approximate search avoiding octree construction
        debugCell_ = mesh_->findCell(debugPoint, polyMesh::FACE_PLANES);

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

        // If reference cell no found use octree search
        // with cell tet-decompositoin
        if (sumHasRef != 1)
        {
            debugCell_ = mesh_->findCell(debugPoint);

            hasRef = (debugCell_ >= 0 ? 1 : 0);
            sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
        }

        if (sumHasRef != 1)
        {
            WarningInFunction
                << "Unable to set debug cell at point " << debugPoint << ":"
                << nl << "Found on " << sumHasRef << " domains (should be one)"
                << endl;
            cellDebugging_ = false;
        }
        else
        {
            cellDebugging_ = true;
        }
    }
    else
    {
        cellDebugging_ = false;
    }
    if (cellDebugging_)
    {
        if (debugCell_ >= 0)
        {
            Pout << "Found debug cell " << debugCell_ << nl << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    addToRunTimeSelectionTable
    (
        solverModule,
        hisaModule,
        dictionary
    );
}

// ************************************************************************* //

