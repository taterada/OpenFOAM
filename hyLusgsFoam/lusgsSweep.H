{
    Field<scalar> D( mesh.V()/dt );
    Field<scalar> D_vib( mesh.V()/dt );

    const volScalarField a( sqrt( gamma*(p/rho)) );
   
    const volScalarField nuMax( "nuMax", lusgsOmegaVis*max(4./3.0 * (mu/rho), alphaTr/rho) );
    
    const volScalarField nuMax_vib( "nuMax_vib", lusgsOmegaVib*max(4./3.0 * (mu/rho), alphaV/rho ) );
    
    Field<scalar> dRhoStar( dRho );
    Field<vector> dRhoUStar( dRhoU );
    Field<scalar> dRhoEStar( dRhoE );
    Field<scalar> dRhoEvStar( dRhoEv );
    
    #include "forwardSweep.H"
    #include "backwardSweep.H"
}
