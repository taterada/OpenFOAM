#include "createMRF.H"

// Primitive variables

Info<< "Reading primitive fields\n" << endl;

volScalarField p // static pressure
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T // translational-rotational temperature
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Tv // vibrational temperature
(
    IOobject
    (
        "Tv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U // velocity
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Thermophysical variables

Info<< "Reading thermophysical properties\n" << endl;

volScalarField R // universal gas constant
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,dimensionedScalar("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), 8.31446261815324)
);

volScalarField M_N2 // N2 molar mass
(
    IOobject
    (
        "M_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,dimensionedScalar("M_N2", dimensionSet(1, 0, 0, 0, -1, 0, 0), 28.014e-3)
);

volScalarField r // specific gas constant
(
    IOobject
    (
        "r",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    R/M_N2
);


volScalarField Cvt // translational specific heat capacity 
(
    IOobject
    (
        "Cvt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.5*r
);

volScalarField Cvr // rotational specific heat capacity
(
    IOobject
    (
        "Cvr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    r
);


volScalarField theta_N2
(
    IOobject
    (
        "theta_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,dimensionedScalar("theta_N2", dimensionSet(0, 0, 0, 1, 0, 0, 0), 3371.)
);

volScalarField ev // vibrational energy
(
    IOobject
    (
        "ev",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (R/M_N2)*theta_N2/(Foam::exp(theta_N2/Tv) - 1.) 
);

volScalarField gamma // specific heat capacity ratio
(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,dimensionedScalar("gamma", dimless, 1.4)
);


//dimensionedScalar As = dimensionedScalar("As", dimensionSet(1, -1, -1, -0.5, 0, 0, 0), 1.401e-06);
//dimensionedScalar Ts = dimensionedScalar("Ts", dimensionSet(0, 0, 0, 1, 0, 0, 0), 107.);

dimensionedScalar A = dimensionedScalar("A", dimless, 2.68e-2);
dimensionedScalar B = dimensionedScalar("B", dimless, 3.18e-1);
dimensionedScalar C = dimensionedScalar("C", dimless, -1.13e1);

dimensionedScalar atm = dimensionedScalar("atm", dimensionSet(1, -1, -2, 0, 0, 0, 0), 101325.);
dimensionedScalar mut = dimensionedScalar("mut", dimless, 14.);

dimensionedScalar avogadro = dimensionedScalar("avogadro", dimless, 6.02214076e23);

volScalarField p_dimless
(
    IOobject
    (
        "p_dimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p/atm
);

dimensionedScalar T_one = dimensionedScalar("T_one", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1.);
dimensionedScalar rho_one = dimensionedScalar("rho_one", dimensionSet(1, -3, 0, 0, 0, 0, 0), 1.);

volScalarField T_dimless
(
    IOobject
    (
        "T_dimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T/T_one
);

dimensionedScalar A_12 = dimensionedScalar("A_12", dimless, 221.53);
dimensionedScalar B_12 = dimensionedScalar("B_12", dimless, 0.029);
dimensionedScalar offset = dimensionedScalar("offset", dimless, 18.42);


dimensionedScalar mu_Blottner = dimensionedScalar("mu_Blottner", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.1);
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu_Blottner*Foam::exp((A*Foam::log(T_dimless) + B)*Foam::log(T_dimless) + C)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu*((5./2.)*Cvt + Cvr)
);
volScalarField alphaTr
(
    IOobject
    (
        "alphaTr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    kappa/(Cvt + Cvr)
);

volScalarField alphaV
(
    IOobject
    (
        "alphaV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.2*mu
);

dimensionedScalar time = dimensionedScalar("time", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.);

volScalarField tau_LT
(
    IOobject
    (
        "tau_LT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    time*(1./p_dimless)*Foam::exp(A_12*(Foam::pow(T_dimless, -1.0/3.0) - B_12) - offset)
);





Info<< "Computing further primitive fields\n" << endl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p/(r*T)
);

volScalarField rho_dimless
(
    IOobject
    (
        "rho_dimless",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho/rho_one
);



volScalarField tau_c
(
    IOobject
    (
        "tau_c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    time*(1./(Foam::sqrt((8.*8.31446261815324*T_dimless)/(3.14159265359*28.014e-3))*(1.0e-21*Foam::pow((50000./T_dimless),2)*Foam::pow(0.01,2))*((avogadro/28.014e-3)*rho_dimless)))
);

volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    tau_LT + tau_c
);


volScalarField etr // translational-rotational energy
(
    IOobject
    (
        "etr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Cvt + Cvr)*T
);


volScalarField e // internal energy
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    etr + ev
);

/*volScalarField Cvv
(
    IOobject
    (
        "Cvv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ev/Tv
);

volScalarField T_effective
(
    IOobject
    (
        "T_effective",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    e/(Cvt + Cvr + Cvv)
);*/



volScalarField Qvt
(
    IOobject
    (
        "Qvt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*(((R/M_N2)*theta_N2/(Foam::exp(theta_N2/T) - 1.)) - ev)/tau
);

// Conservative variables

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

volScalarField rhoEv
(
    IOobject
    (
        "rhoEv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*ev
);

volScalarField heatFlux
(
    IOobject
    (
        "heatFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(- alphaTr*fvc::grad(etr) - alphaV*fvc::grad(ev))
);

// Create numeric flux
numericFlux dbnsFlux
(
    p,
    U,
    T,
    ev,
    Cvt,
    Cvr,
    r
);
