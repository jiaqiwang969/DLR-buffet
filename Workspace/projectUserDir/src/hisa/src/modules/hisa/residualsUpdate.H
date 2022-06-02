/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

    // Construct cont eqn: [A] \Delta W = [R]
    rhoR =
    (
        -fvc::div(phi)
    );

    // Construct mom eqn: [A] \Delta W = [R]
    rhoUR =
    (
        -fvc::div(phiUp)
    );

    // Construct energ eqn: [A] \Delta W = [R]
    rhoER =
    (
        -fvc::div(phiEp)
    );

    // Viscous flow
    if (!inviscid_)
    {
        volScalarField muEff("muEff", turbulence_->muEff());
#if FOUNDATION >= 8
        volScalarField alphaEff("alphaEff", thermophysicalTransport_->alphaEff());
#else
        volScalarField alphaEff("alphaEff", turbulence_->alphaEff());
#endif

        volTensorField& tauMC = tauMC_();
        tauMC = muEff*dev2(Foam::T(fvc::grad(U)));
        rhoUR += fvc::laplacian(muEff,U);
        rhoUR += fvc::div(tauMC);

        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
              fvc::interpolate(muEff)*fvc::snGrad(U)
              + (mesh.Sf()/mesh.magSf() & fvc::interpolate(tauMC))
            )
            & Up
        );
        rhoER += fvc::div(sigmaDotU);

        volScalarField eCalc("eCalc", rhoE/rho-0.5*magSqr(U)); // prevent e inheriting BC from T in thermo
        rhoER += fvc::laplacian(alphaEff, eCalc, "laplacian(alphaEff,e)");
    }

    #include "cellDebug.H"

    // Pseudo and real-time contribution

    tmp<fvScalarMatrix> ddtRho = fvm::ddt(rho);
    tmp<fvVectorMatrix> ddtRhoU = fvm::ddt(rhoU);
    tmp<fvScalarMatrix> ddtRhoE = fvm::ddt(rhoE);
    
#if OPENFOAM >= 1712
    // Work around bug causing timeIndex not to be stored for ddt0 field (CrankNicolson)
    // Seems to be caused by commented out line in GeometricField::storeOldTimes(),
    // commit 44ce956798
    if (mesh.foundObject<volScalarField>("ddt0(rho)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rho)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volVectorField>("ddt0(rhoU)"))
    {
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("ddt0(rhoU)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volScalarField>("ddt0(rhoE)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rhoE)")
        ).timeIndex() = mesh.time().timeIndex();
    }
#endif

    if (!steadyState_)
    {
        // NOTE: Solving the system A \Delta W^\tau = B and therefore the
        // source term needs to be adjusted as below.
        // For steady state analysis this term reduces it to zero.
        //
        // NOTE: The residual is in strong form and therefore the division by V
        rhoR.primitiveFieldRef() -= (ddtRho().diag()*rho-ddtRho().source())/mesh.V();
        rhoUR.primitiveFieldRef() -= (ddtRhoU().diag()*rhoU-ddtRhoU().source())/mesh.V();
        rhoER.primitiveFieldRef() -= (ddtRhoE().diag()*rhoE-ddtRhoE().source())/mesh.V();
    }

    // NOTE: These should be equal if ddtScheme is consistent
    ddtCoeff = max(max(ddtRho().diag(),ddtRhoU().diag()),ddtRhoE().diag())/mesh.V();
