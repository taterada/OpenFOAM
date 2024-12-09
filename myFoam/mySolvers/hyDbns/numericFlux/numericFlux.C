/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

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

#include "fvCFD.H"
#include "numericFlux.H"
#include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::numericFlux::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& ev,
    const volScalarField& Cvt,
    const volScalarField& Cvr,
    const volScalarField& r
)
:
    pFlux_( Foam::dbnsFlux::New( 
        p.mesh(), 
        p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    ),
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    ev_(ev),
    Cvt_(Cvt),
    Cvr_(Cvr),
    r_(r),
    rhoFlux_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate((p_/(r_*T_))*U_) & mesh_.Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate((Cvt_+Cvr_)*T_ + ev_ + 0.5*magSqr(U_))
    ), 
    rhoEvFlux_
    (
        IOobject
        (
            "rhoEvFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(ev_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::numericFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField p_pos( interpolate(p_, pos_) );
    surfaceScalarField p_neg( interpolate(p_, neg_) );

    surfaceVectorField U_pos( interpolate(U_, pos_) );
    surfaceVectorField U_neg( interpolate(U_, neg_) );

    surfaceScalarField T_pos( interpolate(T_, pos_) );
    surfaceScalarField T_neg( interpolate(T_, neg_) );
    
    surfaceScalarField ev_pos( interpolate(ev_, pos_) );
    surfaceScalarField ev_neg( interpolate(ev_, neg_) );

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        // calculate fluxes with reconstructed primitive variables at faces
	pFlux_ -> evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            rhoEvFlux_[faceI],
            p_pos[faceI],  p_neg[faceI],
            U_pos[faceI],  U_neg[faceI],
            T_pos[faceI],  T_neg[faceI],
            ev_pos[faceI], ev_neg[faceI],
            r_[own],       r_[nei],
            Cvt_[own],     Cvt_[nei],
            Cvr_[own],     Cvr_[nei],
            Sf[faceI],
            magSf[faceI]
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pRhoEvFlux = rhoEvFlux_.boundaryFieldRef()[patchi];

        const scalarField& pCvt = Cvt_.boundaryField()[patchi];
        const scalarField& pCvr = Cvr_.boundaryField()[patchi];
        const scalarField& pR  = r_.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pev_pos = ev_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pev_neg = ev_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                pFlux_ -> evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRhoEvFlux[facei],

                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    pT_pos[facei],  pT_neg[facei],
                    pev_pos[facei],  pev_neg[facei],

                    pR[facei],  pR[facei],
                    pCvt[facei], pCvt[facei],
                    pCvr[facei], pCvr[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
        else
        {
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];
            const scalarField& pev = ev_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                // Calculate fluxes
                pFlux_ -> evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRhoEvFlux[facei],
                    
                    pp[facei],  pp[facei],
                    pU[facei],  pU[facei],
                    pT[facei],  pT[facei],
                    pev[facei],  pev[facei],
                    
                    pR[facei],  pR[facei],
                    pCvt[facei], pCvt[facei],
                    pCvr[facei], pCvr[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
    }
}


// ************************************************************************* //
