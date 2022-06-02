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

Application
    getPatchPoints

Description
    Outputs the initial points of the boundary patch
    in the order required for mode shapes
    to the file constant/modal/<patchName>Points.
    e.g. getPatchPoints -patch wing

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Outputs boundary points for patch."
    );

    argList::addOption
    (
        "patch",
        "word",
        "patch name"
    );

   argList::addOption
    (
        "region",
        "word",
        "region to use"
    );

    #include "setRootCase.H"

#if OPENFOAM >= 1906
    if (!args.found("patch"))
#else
    if (!args.optionFound("patch"))
#endif
    {
        FatalErrorIn(args.executable())
            << "Option 'patch' is required."
            << exit(FatalError);
    }

    #include "createTime.H"

    // Create mesh
#if OPENFOAM >= 1906
    const word region = 
        args.getOrDefault<word>("region", Foam::polyMesh::defaultRegion);
#else
    const word region = 
        args.optionLookupOrDefault<word>
        (
            "region", 
            Foam::polyMesh::defaultRegion
        );
#endif

    Foam::Info
    << "Create mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Read initial points from relevant boundary
#if OPENFOAM >= 1906
    word patchName(args.opt<word>("patch"));
#else
    word patchName(args.optionRead<word>("patch"));
#endif
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const vectorField& initialPoints = mesh.boundaryMesh()[patchID].localPoints();

#if OPENFOAM >= 1906
    string fileName = 
        "constant/" + 
        (args.found("region") ? region + "/" : "") + 
        "modal/" + patchName + "Points";
#else
    string fileName = 
        "constant/" + 
        (args.optionFound("region") ? region + "/" : "") + 
        "modal/" + patchName + "Points";
#endif
    OFstream os(fileName);
    IOobject::writeBanner(os) << nl;
    os.writeKeyword("initialPoints") << initialPoints << token::END_STATEMENT << nl;
    IOobject::writeEndDivider(os);

    if (os.good())
    {
        Info << "Successfully wrote points to " << fileName
             << " for patch '" << patchName << "'." << nl << endl;
    }
    else
    {
        FatalErrorIn(args.executable()) << "Could not write points to " << fileName << endl << exit(FatalError);
    }

    return 0;
}
