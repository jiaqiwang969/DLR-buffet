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

Application
    findRoot

Description
    root-finding functions

\*---------------------------------------------------------------------------*/

#include "findRoot.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
rootFunction::rootFunction(double a1, double a2, double a3,double a4)
{ coef1=a1;
  coef2=a2;
  coef3=a3;
  coef3=a4;
}
rootFunction::rootFunction()
{}
    double rootFunction::operator()(const double& x) const // returns the value of the function
    {
        double s1,s2,g1,g2,landa,ret;
        landa=coef2*sqr(x)/coef3;
        s1=1.0+0.275*(1.0-Foam::exp(-35.0*landa))*Foam::exp(-coef1/0.5);
        s2=1.0-(-12.986*landa-123.66*sqr(landa)-405.689*Foam::pow(landa,3.0))*Foam::exp(-Foam::pow(coef1/1.5,1.5));
        ret=coef4/coef3*x;
        if (coef1 <=1.3)
           { g1=1173.51-589.428*coef1+0.2196/(coef1*coef1);
               if (coef2>0.0)
                { return ret-g1*s1;}
             else 
                { return ret-g1*s2;}
            }        
        else
           { g2=331.5*Foam::pow(coef1-0.5658,-0.671);
             if (coef2>0.0)
                { return ret-g2*s1;}
             else 
                { return ret-g2*s2;}
            }
    }

    double rootFunction::d(const double& x) const
    {
        double ds1,ds2,g1,g2,landa,dret,dlanda;
        landa=coef2*sqr(x)/coef3;
        dlanda=2.0*x*coef2/coef3;       
        ds1=0.275*0.35*Foam::exp(-35.0*landa)*Foam::exp(-coef1/0.5)*dlanda;
        ds2=-(-12.986-2.0*123.66*landa-3.0*405.689*Foam::pow(landa,2.0))*Foam::exp(-Foam::pow(coef1/1.5,1.5))*dlanda;
        dret=coef4/coef3;
        if (coef1 <=1.3)
           { g1=1173.51-589.428*coef1+0.2196/(coef1*coef1);
             if (coef2>0.0)
                { return dret-g1*ds1;}
             else 
                { return dret-g1*ds2;}
            }        
        else
           { g2=331.5*Foam::pow(coef1-0.5658,-0.671);
             if (coef2>0.0)
                { return dret-g2*ds1;}
             else 
                { return dret-g2*ds2;}
            }  
    }

    void rootFunction::modify(double x,double y,double z,double w)
    {coef1=max(x,0.027); // Tu
     coef2=y; // duds
     coef3=z; // nu
     coef4=w; // free stream velocity
    }




// ************************************************************************* //
