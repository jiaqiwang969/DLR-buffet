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

Description
    Set the timestep according to the Courant Number / maximum Courant
    number.

    Scale the pseudo Courant number using the Switched Evolution Relaxation
    scheme, by evaluating the change in the L2-norm of the residual (See
    Blazek p 213).

    It is suggested that it increases by a maximum of a factor 2 and decreases
    by not more than a factor of 10.


\*---------------------------------------------------------------------------*/

#include "hisaModule.H"
#include "fvcSmooth.H"


namespace Foam
{

void hisaModule::setPseudoDeltaT()
{
    const fvMesh& mesh = mesh_();
    volScalarField& rPseudoDeltaT = rPseudoDeltaT_();
    {
        tmp< volScalarField > gamma = pThermo_->gamma();
        tmp< volScalarField > psi = pThermo_->psi();
        const volScalarField& ppsi = pThermo_->psi();
        tmp< surfaceScalarField > lambda;
        if (mesh.moving())
        {
            //lambda = fvc::interpolate(sqrt(gamma()/psi())) + mag((fvc::interpolate(U_())&mesh.Sf())-fvc::meshPhi(U_()))/mesh.magSf();
            lambda = fvc::interpolate(sqrt(gamma()/ppsi)) + mag(fvc::interpolate(U_())-fvc::meshPhi(U_())*mesh.Sf()/sqr(mesh.magSf()));
        }
        else
        {
            // lambda = fvc::interpolate(sqrt(gamma()/psi())) + mag(fvc::interpolate(U_())&mesh.Sf()/mesh.magSf());
            lambda = fvc::interpolate(sqrt(gamma()/ppsi)) + fvc::interpolate(mag(U_()));
        }

        if (localTimestepping_)
        {
            rPseudoDeltaT == dimensionedScalar("0", dimless/dimTime, 0.0);
            // surfaceScalarField frdt = mesh.deltaCoeffs()*lambda;
            surfaceScalarField frdt = mesh.nonOrthDeltaCoeffs()*(lambda());
            const labelList& own = mesh.owner();
            const labelList& nei = mesh.neighbour();
            forAll(own,facei)
            {
                rPseudoDeltaT[own[facei]] =
                    max(rPseudoDeltaT[own[facei]],frdt[facei]);
                rPseudoDeltaT[nei[facei]] =
                    max(rPseudoDeltaT[nei[facei]],frdt[facei]);
            }
            forAll(frdt.boundaryField(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                if (mesh.boundary()[patchi].coupled())
                {
                    forAll(fc, bfacei)
                    {
                        rPseudoDeltaT[fc[bfacei]] =
                            max
                            (
                                rPseudoDeltaT[fc[bfacei]],
                                frdt.boundaryField()[patchi][bfacei]
                            );
                    }
                }
                else if 
                (
                    mesh.boundary()[patchi].type() == "wall"
                )
                {
                    const scalarField pNonOrthDeltaCoef(mesh.nonOrthDeltaCoeffs().boundaryField()[patchi]);

                    forAll(fc, bfacei)
                    {
                        // Value based on cell value lambda and twice the patch edge length
                        scalar pLambda = 
                            0.5*pNonOrthDeltaCoef[bfacei]*
                            (
                                sqrt(gamma()[fc[bfacei]]/psi()[fc[bfacei]]) + 
                                mag(U_()[fc[bfacei]])
                            );
                    
                        rPseudoDeltaT[fc[bfacei]] =
                            max
                            (
                                rPseudoDeltaT[fc[bfacei]], // Cell value based on max of internal neighbours
                                pLambda 
                            );
                    }
                }
            }
            gamma.clear();
            psi.clear();
            lambda.clear();

            rPseudoDeltaT /= pseudoCoField_();
            rPseudoDeltaT.correctBoundaryConditions();
        }
        else
        {
            rPseudoDeltaT = max(mesh.deltaCoeffs()*lambda); // Global timestep
            rPseudoDeltaT /= pseudoCoNum_();
        }
    }

    if (localTimestepping_)
    {
        scalar totCells = mesh.globalData().nTotalCells();
        Info<< "Pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField_()).value()
            << " Mean: "
            << sum(pseudoCoField_()).value()/totCells
            << " Max: "
            << max(pseudoCoField_()).value() << endl;
        Info<< "Pseudo deltaT: "
            << "Min: "
            << 1.0/max(rPseudoDeltaT).value()
            << " Mean: "
            << sum(1.0/rPseudoDeltaT).value()/totCells
            << " Max: "
            << 1.0/min(rPseudoDeltaT).value() << endl;
    }
    else
    {
        Info << "Pseudo Courant No: " << pseudoCoNum_->value() << endl;
        Info << "Pseudo deltaT: " << 1.0/min(rPseudoDeltaT).value() << endl;
    }
}


void hisaModule::setPseudoCoNum()
{
    // Switched Evolution Relaxation
    if (initRes_.valid())
    {
        if (!solnControl_->firstIter() && prevRes_.valid())
        {
            residualIO& initRes(initRes_());
            residualIO& prevRes(prevRes_());

            scalar normInit = sqrt(magSqr(initRes.getScalar(0)) +
                                magSqr(initRes.getScalar(1)) +
                                magSqr(initRes.getVector(0)));
            scalar normPrev = sqrt(magSqr(prevRes.getScalar(0)) +
                                magSqr(prevRes.getScalar(1)) +
                                magSqr(prevRes.getVector(0)));
            scalar coNumRatio = normPrev/normInit;
            coNumRatio = max(min(coNumRatio, pseudoCoNumMaxIncr_), pseudoCoNumMinDecr_);
            if (localTimestepping_)
            {
                pseudoCoField_() *= coNumRatio;
                pseudoCoField_() = max(min(pseudoCoField_(), pseudoCoNumMax_), pseudoCoNumMin_);
                pseudoCoField_().correctBoundaryConditions();
            }
            else
            {
                pseudoCoNum_() *= coNumRatio;
                pseudoCoNum_().value() = max(min(pseudoCoNum_().value(), pseudoCoNumMax_), pseudoCoNumMin_);
            }
        }
        prevRes_.clear();
        prevRes_.set(new residualIO(initRes_()));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
