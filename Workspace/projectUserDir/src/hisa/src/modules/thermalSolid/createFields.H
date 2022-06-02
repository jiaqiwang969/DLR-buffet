// Populate solid field pointer lists

Info<< "Creating solid thermo\n" << endl;
thermo_ = solidThermo::New(mesh);
solidThermo& thermo = thermo_();

//Info<< "Creating radiation model\n" << endl;
//radiation_.set(i, radiation::radiationModel::New(thermo.T()));

Info<< "Reading fvOptions\n" << endl;
solidHeatSources_.set
(
#if FOUNDATION >= 9
    new fvModels(mesh)
#else
    new fv::options(mesh)
#endif
);
#if FOUNDATION >= 9
solidConstraints_.set
(
    new fvConstraints(mesh)
);
#endif

if (!thermo.isotropic())
{
    Info<< "Adding coordinateSystem\n" << endl;
#if FOUNDATION >= 9
    coordinates_ = coordinateSystem::New(mesh, thermo.properties());
#else
    coordinates_ = coordinateSystem::New(mesh, thermo);
#endif    

    tmp<volVectorField> tkappaByCp =
        thermo.Kappa()/thermo.Cp();

    aniAlpha_.set
    (
        new volSymmTensorField
        (
            IOobject
            (
                "Anialpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            (
                "zero",
                tkappaByCp().dimensions(),
                Zero
            )
        )
    );

    updateAniAlpha(tkappaByCp);
}

IOobject betavSolidIO
(
    "betavSolid",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

#ifdef FOAM_VERSION_4
if (betavSolidIO.headerOk())
#else
if (betavSolidIO.typeHeaderOk<volScalarField>())
#endif
{
    betavSolid_.set
    (
        new volScalarField(betavSolidIO, mesh)
    );
}
else
{
    betavSolid_.set
    (
        new volScalarField
        (
            IOobject
            (
                "betavSolid",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("1", dimless, scalar(1.0))
        )
    );
}
