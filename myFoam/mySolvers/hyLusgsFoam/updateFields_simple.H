{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

    // Calculate internal energy from rhoE

    e = (rhoE/rho - 0.5*magSqr(U));
    
    // Calculate vibrational energy from rhoEv
    ev = rhoEv/rho;
    
    dimensionedScalar evMin = 0.5*(min(R/M_N2))*min(theta_N2)/(Foam::exp(min(theta_N2)/TvMin) - 1.);
    dimensionedScalar evMax = (max(R/M_N2))*max(theta_N2)/(Foam::exp(max(theta_N2)/TvMax) - 1.);
    
    // Bound vibrational energy
    boundMinMax(ev, evMin, evMax);
    
    dimensionedScalar CvMin = min(Cvt+Cvr);
    dimensionedScalar CvMax = max(Cvt+Cvr);

    dimensionedScalar eMin = CvMin*TMin + min(ev);
    dimensionedScalar eMax = CvMax*TMax + max(ev);
	
    // Bound internal energy
    boundMinMax(e, eMin, eMax);
    
    etr = e - ev;

    Tv = theta_N2/Foam::log(r*theta_N2/ev + 1.);
    Tv.correctBoundaryConditions();
    
    // Compute T from e and ev
    T = (etr)/(Cvt + Cvr);
    T.correctBoundaryConditions();
    
    // Compute p from rho and T
    p = rho*r*T;
    p.correctBoundaryConditions();
    
    // Bound pressure
    boundMinMax(p, pMin, pMax);
   
    
    // Update boundary field of rho
    rho.boundaryFieldRef() = p.boundaryField() / ( r.boundaryField() * T.boundaryField() );
    // Bound density
    boundMinMax(rho, rhoMin, rhoMax);
    
    ev.boundaryFieldRef() = (R.boundaryField()/M_N2.boundaryField())*theta_N2.boundaryField()/(Foam::exp(theta_N2.boundaryField()/Tv.boundaryField()) - 1.); 
    e.boundaryFieldRef() = (Cvt.boundaryField()+Cvr.boundaryField())*T.boundaryField() + ev.boundaryField();
    
    etr.boundaryFieldRef() = e.boundaryField() - ev.boundaryField();
     
    // Update boundary field of rhoU 
    rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() = rho.boundaryField()*( e.boundaryField() + 0.5*magSqr(U.boundaryField()) );
    
    // Update boundary field of rhoEv
    rhoEv.boundaryFieldRef() = rho.boundaryField()*ev.boundaryField();
    
}
