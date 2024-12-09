{
    //Cvv = ev/Tv;
    //Cvv.boundaryFieldRef() = ev.boundaryField()/Tv.boundaryField();
    
    //T_effective = e/(Cvt + Cvr + Cvv);
    //T_effective.boundaryFieldRef() = e.boundaryField()/(Cvt.boundaryField() + Cvr.boundaryField() + Cvv.boundaryField());
    
    p_dimless = p/atm;
    p_dimless.boundaryFieldRef() = p.boundaryField() / atm.value();
    
    
    T_dimless = T/T_one;
    T_dimless.boundaryFieldRef() = T.boundaryField()/T_one.value();
    
    rho_dimless = rho/rho_one;
    rho_dimless.boundaryFieldRef() = rho.boundaryField()/rho_one.value();
    
    mu = mu_Blottner*Foam::exp((A*Foam::log(T_dimless) + B)*Foam::log(T_dimless) + C);
    mu.boundaryFieldRef() = mu_Blottner.value()*Foam::exp((A.value()*Foam::log(T_dimless.boundaryField()) + B.value())*Foam::log(T_dimless.boundaryField()) + C.value());
    
    kappa = mu*((5./2.)*Cvt + Cvr);
    kappa.boundaryFieldRef() = mu.boundaryField()*((5./2.)*Cvt.boundaryField() + Cvr.boundaryField());
    
    alphaTr = kappa/(Cvt + Cvr);
    alphaTr.boundaryFieldRef() = kappa.boundaryField()/(Cvt.boundaryField() + Cvr.boundaryField());
    
    alphaV = 1.2*mu;
    alphaV.boundaryFieldRef() = 1.2*mu.boundaryField();
       
       
    //TODO: tadu je to s tou efektivni teplotou celkem na vazkach
    tau_LT = time*(1./p_dimless)*Foam::exp(A_12*(Foam::pow(T_dimless, -1.0/3.0) - B_12) - offset);
    tau_LT.boundaryFieldRef() = time.value()*(1./p_dimless.boundaryField())*Foam::exp(A_12.value()*(Foam::pow(T_dimless.boundaryField(), -1.0/3.0) - B_12.value()) - offset.value());
    
    tau_c = time*(1./(Foam::sqrt((8.*8.31446261815324*T_dimless)/(3.14159265359*28.014e-3))*(1.0e-21*Foam::pow((50000./T_dimless),2)*Foam::pow(0.01,2))*((avogadro/28.014e-3)*rho_dimless)));
    tau_c.boundaryFieldRef() = time.value()*(1./(Foam::sqrt((8.*8.31446261815324*T_dimless.boundaryField())/(3.14159265359*28.014e-3))*(1.0e-21*Foam::pow((50000./T_dimless.boundaryField()),2)*Foam::pow(0.01,2))*((avogadro.value()/28.014e-3)*rho_dimless.boundaryField())));
    
    tau = tau_LT + tau_c;
    tau.boundaryFieldRef() = tau_LT.boundaryField() + tau_c.boundaryField();

    Qvt = rho*(((R/M_N2)*theta_N2/(Foam::exp(theta_N2/T) - 1.)) - ev)/tau;
    Qvt.boundaryFieldRef() = rho.boundaryField()*(((R.boundaryField()/M_N2.boundaryField())*theta_N2.boundaryField()/(Foam::exp(theta_N2.boundaryField()/T.boundaryField()) - 1.)) - ev.boundaryField())/tau.boundaryField();
    
    heatFlux = mag(- alphaTr*fvc::grad(etr)- alphaV*fvc::grad(ev));
    
    
}
