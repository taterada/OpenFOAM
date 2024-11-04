/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
    
\*---------------------------------------------------------------------------*/

#include "hllLMFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(hllLMFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, hllLMFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllLMFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    scalar& rhoEvFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& evLeft,
    const scalar& evRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvtLeft,
    const scalar& CvtRight,
    const scalar& CvrLeft,
    const scalar& CvrRight,
    const vector& Sf,
    const scalar& magSf
) const
{
  //if (mag(meshPhi)>0.0) 
  //  {
  //    FatalError
  //      << "This dbnsFlux is not ready to run with moving meshes." << nl
  //      << exit(FatalError);
  //  };

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = 7./5.;
    const scalar kappaRight = 7./5.;

    // Compute conservative variables assuming perfect gas law

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*((CvtLeft+CvrLeft)*TLeft+evLeft+0.5*magSqr(ULeft));
    const scalar rhoERight = rhoRight*((CvtRight+CvrRight)*TRight+evRight+0.5*magSqr(URight));
    
    const scalar rhoEvLeft = rhoLeft*evLeft;
    const scalar rhoEvRight = rhoRight*evRight;

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector);
    const scalar qRight = (URight & normalVector);

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(0.0,kappaLeft * pLeft/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0.0,kappaRight * pRight/rhoRight));

    // Compute Roe weights
    const scalar rhoLeftSqrt = Foam::sqrt(max(0.0,rhoLeft));
    const scalar rhoRightSqrt = Foam::sqrt(max(0.0,rhoRight));

    const scalar wLeft = rhoLeftSqrt
        /stabilise((rhoLeftSqrt + rhoRightSqrt),VSMALL);

    const scalar wRight = 1 - wLeft;

    // Roe averaged velocity
    const vector UTilde = wLeft*ULeft + wRight*URight;

    // Roe averaged contravariant velocity
    const scalar contrUTilde = (UTilde & normalVector);

    // Roe averaged total enthalpy
    const scalar HTilde = wLeft*HLeft + wRight*HRight;

    // Roe averaged kappa
    const scalar kappaTilde = wLeft*kappaLeft + wRight*kappaRight;

    // Speed of sound with Roe reconstruction values
    const scalar aTilde =
        Foam::sqrt(max(0 ,(kappaTilde - 1)*(HTilde - 0.5*sqr(contrUTilde))));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft - aLeft, contrUTilde - aTilde);
    const scalar SRight = max(contrUTilde + aTilde, qRight + aRight);



    if (pos(SLeft))
    {
        rhoFlux  = (qLeft*rhoLeft)*magSf;
        rhoUFlux = (qLeft*rhoULeft + pLeft*normalVector)*magSf;
        rhoEFlux = (qLeft*(rhoELeft + pLeft))*magSf;
        rhoEvFlux  = (qLeft*rhoEvLeft)*magSf;
    }
    else if (neg(SRight))
    {
        rhoFlux  = (qRight*rhoRight)*magSf;
        rhoUFlux = (qRight*rhoURight + pRight*normalVector)*magSf;
        rhoEFlux = (qRight*(rhoERight + pRight))*magSf;
        rhoEvFlux  = (qLeft*rhoEvRight)*magSf;
    }
    else
    {
        rhoFlux  = (SRight*(qLeft*rhoLeft) - SLeft*(qRight*rhoRight))*magSf/(SRight - SLeft);
        rhoUFlux = (SRight*(qLeft*rhoULeft + pLeft*normalVector) -
                    SLeft*(qRight*rhoURight + pRight*normalVector))*magSf/(SRight - SLeft);
        rhoEFlux = (SRight*(qLeft*(rhoELeft + pLeft)) -
                    SLeft*(qRight*(rhoERight + pRight)))*magSf/(SRight - SLeft);
        rhoEvFlux  = (SRight*(qLeft*rhoEvLeft) - SLeft*(qRight*rhoEvRight))*magSf/(SRight - SLeft);

        // Compute Mach number function and adjusted pressure in star region
        const scalar MachLeft = mag(ULeft)/aLeft;
        const scalar MachRight = mag(URight)/aRight;
        const scalar LIM = min(max(MachLeft,MachRight),1.0);

        rhoFlux += SRight*SLeft/(SRight - SLeft)*(rhoRight - rhoLeft)*magSf;
        rhoUFlux += LIM*SRight*SLeft/(SRight - SLeft)*(rhoURight - rhoULeft)*magSf;
        rhoEFlux += SRight*SLeft/(SRight - SLeft)*(rhoERight - rhoELeft)*magSf;
        rhoEvFlux += SRight*SLeft/(SRight - SLeft)*(rhoEvRight - rhoEvLeft)*magSf;
    }

}


// ************************************************************************* //
