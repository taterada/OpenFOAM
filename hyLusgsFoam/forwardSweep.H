{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];
            scalar ac = 0.5 * lusgsOmegaInv * ( mag((U[ci] & pSf[facei]))  
					+ a[ci]*mag(pSf[facei]));
            D[ci] += ac;
            D_vib[ci] += ac;
        }
    }

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = mesh.V()[cellI]/dt.value() * dRho[cellI];
        vector rRhoU = mesh.V()[cellI]/dt.value() * dRhoU[cellI];
        scalar rRhoE = mesh.V()[cellI]/dt.value() * dRhoE[cellI];
        scalar rRhoEv = mesh.V()[cellI]/dt.value() * dRhoEv[cellI];


        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index

            if (mesh.isInternalFace(faceI))
            {
                // Calculate D
                label own = owner[faceI];
                label nei = neighbour[faceI];
		        scalar dvol = mag( (mesh.C()[own] - mesh.C()[nei]) & Sf[faceI] );

                scalar ac = 0.5 * lusgsOmegaInv * ( 
					   mag( (U[cellI] & Sf[faceI])) 
					   + a[cellI] * magSf[faceI] );
                scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];
                scalar av_vib = sqr(magSf[faceI]) / dvol * nuMax_vib[cellI];

                D[cellI] += ac + av;
                D_vib[cellI] += ac + av_vib;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    scalar ac = 0.5 * lusgsOmegaInv * (mag( (U[own] & Sf[faceI])) + a[own] * magSf[faceI] );
                    scalar av = sqr(magSf[faceI]) / dvol * nuMax[own];
                    scalar av_vib = sqr(magSf[faceI]) / dvol * nuMax_vib[own];
                    
                    scalar rho1  = rho[own] + dRhoStar[own];
                    vector rhoU1 = rhoU[own] + dRhoUStar[own];
                    scalar rhoE1 = rhoE[own] + dRhoEStar[own];
                    scalar rhoEv1 = rhoEv[own] + dRhoEvStar[own];
                    
                    scalar p1    = (r[own]/(Cvt[own]+Cvr[own])) * (rhoE1 - rhoEv1 - 0.5*magSqr(rhoU1)/rho1);
		    
		    
                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];

                    rRho  += (ac + av) * dRhoStar[own]  + 
                        0.5 * (rho1*phi1 - rho[own]*phi0); 
                    rRhoU += (ac + av) * dRhoUStar[own] +
                        0.5 * ( (rhoU1*phi1 - rhoU[own]*phi0) + (p1 - p[own])*Sf[faceI] );
                    rRhoE += (ac + av) * dRhoEStar[own] +
		                0.5 * (rhoE1*phi1+p1*phi1 - (rhoE[own]*phi0+p[own]*phi0));
		            rRhoEv += (ac + av_vib) * dRhoEvStar[own] +
		                0.5 * (rhoEv1*phi1 - rhoEv[own]*phi0); 
                }
            }
        }

        dRhoStar[cellI]  = rRho / D[cellI];
        dRhoUStar[cellI] = rRhoU / D[cellI];
        dRhoEStar[cellI] = rRhoE / D[cellI];
        dRhoEvStar[cellI] = rRhoEv / D_vib[cellI];
    }
}
