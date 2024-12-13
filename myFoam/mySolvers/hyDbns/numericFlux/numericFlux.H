/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM:
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of OpenFOAM nor foam-extend.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    numericFlux

Description
    Single level numeric flux class for density-based solvers. Based on
    numericFlux by Jemcov and Jasak (foam-extend)

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak
    Rewrite by Jiri Furst

SourceFiles
    numericFlux.H
    numericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef numericFlux_H
#define numericFlux_H

#include "numericFluxBase.H"
#include "dbnsFlux.H"
#include "basicThermo.H"
#include "MRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class numericFlux Declaration
\*---------------------------------------------------------------------------*/

class numericFlux
:
    public numericFluxBase
{
    // Private data

    autoPtr<dbnsFlux> pFlux_;
    
    //- Reference to mesh
    const fvMesh& mesh_;
    
    // Reference to primitive fields
    
    //- Static pressure
    const volScalarField& p_;
    
    //- Velocity
    const volVectorField& U_;
    
    //- Static temperature
    const volScalarField& T_;
    
    const volScalarField& ev_;
    
    const volScalarField& Cvt_;
    const volScalarField& Cvr_;
    const volScalarField& r_;


    // Fluxes
    
    //- Density flux
    surfaceScalarField rhoFlux_;
    
    //- Velocity flux
    surfaceVectorField rhoUFlux_;
    
    //- Energy flux
    surfaceScalarField rhoEFlux_;
    
    surfaceScalarField rhoEvFlux_;
    
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    numericFlux(const numericFlux&);
    
    //- Disallow default bitwise assignment
    void operator=(const numericFlux&);
    
    
    //- Return internal field of mass flux
    const scalarField& rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }
    
    //- Return access to internal field of mass flux
    //scalarField& rhoFluxI()
    //    {
    //        return rhoFlux_.internalField();
    //    }
    
    //- Return internal field of momentum flux
    const vectorField& rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }
    
    //- Return access to internal field of momentum flux
    //vectorField& rhoUFluxI()
    //    {
    //        return rhoUFlux_.internalField();
    //    }
    
    //- Return access to internal field of energy flux
    const scalarField& rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }
    
    const scalarField& rhoEvFluxI() const
        {
            return rhoEvFlux_.internalField();
        }
    
    //- Return access to internal field of energy flux
    //scalarField& rhoEFluxI()
    //    {
    //        return rhoEFlux_.internalField();
    //    }
    
    
public:
    
    // Constructors
    
    //- Construct from components
    numericFlux
    (
        const volScalarField& p,
        const volVectorField& U,
        const volScalarField& T,
        const volScalarField& ev,
        const volScalarField& Cvt,
        const volScalarField& Cvr,
        const volScalarField& r
    );
    
    
    //- Destructor
    virtual ~numericFlux()
        {}

    
    // Member Functions
    
    //- Return mesh reference
    const fvMesh& mesh() const
        {
            return mesh_;
        }

    
    // Return fluxes
    
    //- Return density flux
    virtual const surfaceScalarField& rhoFlux() const
        {
            return rhoFlux_;
        }
    
    //- Return velocity flux
    virtual const surfaceVectorField& rhoUFlux() const
        {
            return rhoUFlux_;
        }
    
    //- Return energy flux
    virtual const surfaceScalarField& rhoEFlux() const
        {
            return rhoEFlux_;
        }
        
    virtual const surfaceScalarField& rhoEvFlux() const
        {
            return rhoEvFlux_;
        }
    
    
    // Return residuals
    
    //- Return density equation residual
    virtual tmp<scalarField> rhoResidual() const
        {
            return fvc::div(rhoFlux_)().internalField();
        }
    
    //- Return momentum equation flux
    virtual tmp<vectorField> rhoUResidual() const
        {
            return fvc::div(rhoUFlux_)().internalField();
        }
    
    //- Return energy equation flux
    virtual tmp<scalarField> rhoEResidual() const
        {
            return fvc::div(rhoEFlux_)().internalField();
        }
        
    virtual tmp<scalarField> rhoEvResidual() const
        {
            return fvc::div(rhoEvFlux_)().internalField();
        }
        
    // Update fluxes based on current state
    
    //- Compute flux
    virtual void computeFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
