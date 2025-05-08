
    //- get the currernt time and time step
    scalar t = runTime.time().value();
    dimensionedScalar deltaT = runTime.deltaT();


    volScalarField& rhoEff = rhoEffMetalHydride[i];
    volScalarField& CpEff = CpEffMetalHydride[i];
    volScalarField& KappaEff = KappaEffMetalHydride[i];
    volScalarField& pEq = pEqMetalHydride[i];
    volScalarField& SmMH = SmMHMetalHydride[i];
    volScalarField& SeMH = SeMHMetalHydride[i];

    dimensionedScalar& deltaH = deltaHMetalHydride[i];
    dimensionedScalar& deltaS = deltaSMetalHydride[i];
    dimensionedScalar& Ca = CaMetalHydride[i];
    dimensionedScalar& Ea = EaMetalHydride[i];
    dimensionedScalar& rhos = rhosMetalHydride[i];
    dimensionedScalar& pRef = pRefMetalHydride[i];
    dimensionedScalar& R = RMetalHydride[i];

    dimensionedScalar& rhog = rhogMetalHydride[i];
    dimensionedScalar& Cpg = CpgMetalHydride[i];
    dimensionedScalar& kappag = kappagMetalHydride[i];
    dimensionedScalar& Mg = MgMetalHydride[i];

    scalar& a = aMetalHydride[i];
    scalar& b = bMetalHydride[i];
    dimensionedScalar Pgas("Pgas", dimPressure, a*t + b);


    porousZones& pZones = pZonesSolid[i];
    volScalarField& porosity = porositySolid[i];

    
    pEq = pRef*exp((deltaH/(R*T)) - (deltaS/R));
    SmMH = Ca*exp(-1.0*(Ea/(R*T)))*log(Pgas/pEq)*(rhos - rho);

    fvScalarMatrix rhoEqn
    (
        (1.0 - porosity)*fvm::ddt(rho)
      ==
        SmMH
    );

    rhoEqn.relax();
    rhoEqn.solve();


    SeMH = SmMH*(deltaH/Mg);
    rhoEff = rhog*porosity + rho*(1.0 - porosity);
    CpEff = (rhog*Cpg*porosity + rho*cp*(1.0 - porosity))/rhoEff;
    KappaEff = kappag*porosity + Kappa*(1.0 - porosity);


{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rhoEff*CpEff, T)
          - fvm::laplacian(KappaEff, T)
          ==
            SeMH*1e-1
        );

        TEqn().relax();
        TEqn().solve();
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}

