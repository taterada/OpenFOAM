dictionary lusgsDict = mesh.solutionDict().subDict("LUSGS");

label lusgsIntIters = 
    lusgsDict.lookupOrDefault<label>("internalIterations", 1);

scalar lusgsTolerance =
    lusgsDict.lookupOrDefault<scalar>("tolerance", 0.0);

scalar lusgsRelTol =
    lusgsDict.lookupOrDefault<scalar>("relTol", 0.0);

// relaxation factor 1<omega_inv<=2
scalar lusgsOmegaInv = 
     lusgsDict.lookupOrDefault<scalar>("omega_inv", 1.2);
     
// adjustment of the viscous term estimation 1<=omega_vis
scalar lusgsOmegaVis = 
     lusgsDict.lookupOrDefault<scalar>("omega_vis", 1.);
     
scalar lusgsOmegaVib = 
     lusgsDict.lookupOrDefault<scalar>("omega_vib", 1.);

bool LTS = 
    lusgsDict.lookupOrDefault("localTimestepping", false);
