/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at you
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    myLusgsFoam

Description
    Density-based compressible implicit steady-state & transient flow solver
    using LU-SGS method.
    
    \verbatim
        Furst, J. (2018).
        Development of a coupled matrix-free LU-SGS solver for turbulent compressible flows,
        Computers & Fluids, https://doi.org/10.1016/j.compfluid.2018.04.020
    \endverbatim


Author
    Jiri Furst

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "boundMinMax.H"
#include "numericFlux.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readLUSGSControls.H"
    #include "createRDeltaTau.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Patch for correct calculation of meshPhi (U needs old values)
    {
        auto dummy = fvc::ddt(U);
    } 
  
    Info<< "\nStarting time loop\n" << endl;

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
   
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readFieldBounds.H"

        surfaceScalarField amaxSf("amaxSf", 
        mag(fvc::interpolate(U) & mesh.Sf()) +
        mesh.magSf() * fvc::interpolate(sqrt(gamma*(p/rho))));
        

        #include "compressibleCFLNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        Info<< "\n Time = " << runTime.value() << nl;
        
        scalar initialRezRho=0, initialRezRhoU=0, initialRezRhoE=0, initialRezRhoEv=0;

        for (int intIter=0; intIter<lusgsIntIters; intIter++)
        {
            Info << "LUSGS: iteration " << intIter+1 << nl;
            
            // Solve the approximate Riemann problem for this time step
            dbnsFlux.computeFlux();
            
            dimensionedScalar dt = runTime.deltaT();

            volSymmTensorField devRhoReff("devRhoReff", -mu*dev(twoSymm(fvc::grad(U))));

            volScalarField dRho( -dt * ( 
                fvc::ddt(rho) + fvc::div(dbnsFlux.rhoFlux()) 
            ));

            volVectorField dRhoU( -dt * (
                fvc::ddt(rhoU)
                + fvc::div(dbnsFlux.rhoUFlux())
                + fvc::div(devRhoReff) 
            )); 

            volScalarField dRhoE(-dt*( 
                fvc::ddt(rhoE)
                + fvc::div(dbnsFlux.rhoEFlux()) 
                + fvc::div(devRhoReff & U)
                - fvc::laplacian(alphaTr, etr)
                - fvc::laplacian(alphaV, ev) 
            ));

            volScalarField dRhoEv(-dt*( 
                fvc::ddt(rhoEv)
                + fvc::div(dbnsFlux.rhoEvFlux()) 
                - fvc::laplacian(alphaV, ev)
                - Qvt 
            ));

            scalar rezRho  = fvc::domainIntegrate( mag(dRho) / dt ).value();
            scalar rezRhoU = fvc::domainIntegrate( mag(dRhoU) / dt ).value();
            scalar rezRhoE = fvc::domainIntegrate( mag(dRhoE) / dt ).value();
            scalar rezRhoEv = fvc::domainIntegrate( mag(dRhoEv) / dt ).value();

            if (intIter == 0) {
                initialRezRho  = rezRho;
                initialRezRhoU = rezRhoU;
                initialRezRhoE = rezRhoE;
                initialRezRhoEv = rezRhoEv;
            }
            
            #include "lusgsSweep.H"
 
            rho  += dRho;
            rhoU += dRhoU;
            rhoE += dRhoE;
            rhoEv += dRhoEv;
            

            #include "updateFields_simple.H"
            #include "updateThermoFields_simple.H"

	        scalar finalRezRho  = fvc::domainIntegrate( mag(dRho) / dt ).value();
	        scalar finalRezRhoU = fvc::domainIntegrate( mag(dRhoU) / dt ).value();
	        scalar finalRezRhoE = fvc::domainIntegrate( mag(dRhoE) / dt ).value();
	        scalar finalRezRhoEv = fvc::domainIntegrate( mag(dRhoEv) / dt ).value();
	        
	        Info << "LUSGS:  Solving for rho,  " << "Initial residual = " << rezRho << ", " << "Final residual = " << finalRezRho << ", No Iterations 1" << nl;
	        Info << "LUSGS:  Solving for rhoU, " << "Initial residual = " << rezRhoU << ", " << "Final residual = " << finalRezRhoU << ", No Iterations 1" << nl;
	        Info << "LUSGS:  Solving for rhoE, " << "Initial residual = " << rezRhoE << ", " << "Final residual = " << finalRezRhoE << ", No Iterations 1" << nl;
	        Info << "LUSGS:  Solving for rhoEv, " << "Initial residual = " << rezRhoEv << ", " << "Final residual = " << finalRezRhoEv << ", No Iterations 1" << nl;

            bool lastIteration = ( intIter + 1 == lusgsIntIters );

            lastIteration = lastIteration || ( 
                (finalRezRho < lusgsRelTol * initialRezRho ) && 
                (finalRezRhoU < lusgsRelTol * initialRezRhoU ) && 
                (finalRezRhoE < lusgsRelTol * initialRezRhoE ) && 
                (finalRezRhoEv < lusgsRelTol * initialRezRhoEv ) );

            lastIteration = lastIteration || ( 
                (finalRezRho < lusgsTolerance) && 
                (finalRezRhoU < lusgsTolerance) && 
                (finalRezRhoE < lusgsTolerance) && 
                (finalRezRhoEv < lusgsTolerance) );


            if (lastIteration) 
            {
                Info << "LUSGS: converged in " << intIter+1 << " iterations." << nl;
                break;
            }
            
        }
	    runTime.write();
	
        Info<< "    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
	
     
    }

    Info<< "\n end \n";
    
    return(0);
}


// ************************************************************************* //
