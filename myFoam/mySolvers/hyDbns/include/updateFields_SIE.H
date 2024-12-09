{
    // Estimate the bdensity at the boundary using "zeroGradient"
    rho.boundaryFieldRef() = thermo.rho()().boundaryField();

    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();
    MRF.correctBoundaryVelocity(U);

    // Calculate internal energy from rhoE
    const volScalarField Cv = thermo.Cv();

    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();

    // Bound enthalpy
    dimensionedScalar CvMin = min(Cv);
    dimensionedScalar CvMax = max(Cv);

    dimensionedScalar eMin = CvMin*TMin;
    dimensionedScalar eMax = CvMax*TMax;

    boundMinMax(e, eMin, eMax);

    // Compute p from rho
    p = (thermo.Cp()/thermo.Cv() - 1.0) * rho * e;
    p.correctBoundaryConditions(); 

    // Bound pressure
    boundMinMax(p, pMin, pMax);

    // Correct thermodynamics
    thermo.correct();

    // Update boundary field of rho
    rho.boundaryFieldRef() = thermo.psi().boundaryField() * p.boundaryField();

    // Bound density
    boundMinMax(rho, rhoMin, rhoMax);

    // Update boundary field of rhoU 
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() =
      rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()) );
    
}
