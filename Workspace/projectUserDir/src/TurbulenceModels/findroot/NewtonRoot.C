/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    NewtonRoot

Description
    Newton root Based on Numerical Recipes in C++, Section 9.1, page 358.
    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)


\*----------------------------------------------------------------------------*/

#include "NewtonRoot.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Func>
const Foam::label Foam::NewtonRoot<Func>::maxIter = 60;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Func>
Foam::NewtonRoot<Func>::NewtonRoot(const Func& f, const scalar eps)
:
    f_(f),
    eps_(eps)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Func>
Foam::scalar Foam::NewtonRoot<Func>::root
(
    const scalar x0

) const
{
    scalar f, df, resid=1, Iter=1,xstart=x0;
    
    while(Iter<maxIter)
    {
     f = f_(xstart);
     df = f_.d(xstart);
     resid=f_(xstart);
     xstart=xstart-f/df;
     if (mag(resid)<eps_)
        { return xstart;}

     Iter=Iter+1;

    } 
    FatalErrorIn
    (
        "Foam::scalar Foam::NewtonRoot<Func>::root\n"
        "(\n"
        "    const scalar x0,\n"
        "    const scalar x1\n"
        ") const"
    )   << "Maximum number of iterations exceeded" << abort(FatalError);

    // Dummy return to keep compiler happy
    return xstart;
}




// ************************************************************************* //
