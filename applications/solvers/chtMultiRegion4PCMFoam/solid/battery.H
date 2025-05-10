    volScalarField& SeB = SeBattery[i];
    word& state = stateBattery[i];
    scalar& C = CBattery[i];
    dimensionedScalar& R = RBattery[i];
    dimensionedScalar& I = IBattery[i];
    dimensionedScalar& dEdT = dEdTBattery[i];

    //- determine the sign of currernt
    scalar sign;
    if(state == "discharge"){sign = -1.0;}
    else{sign = 1.0;}


    //- get the currernt time and time step
    scalar t = runTime.time().value();
    scalar deltat = runTime.deltaT().value();

    if(t > tOld)
    {
        //- compute the SOC
        scalar SOCCurr = SOC[i] + sign*(100/(3600*C))*I.value()*deltat;
        SOC[i] = SOCCurr;
        if( min(SOC) <= 0.0 )
        {
            Info << "All batteries are discharged" << endl;
            FatalErrorIn("finished") << exit(FatalError);
        }
        Info << "SOC: " << SOC[i] << endl;
        tOld = t;
    }

    dimensionedScalar V = sum(solidRegions[i].V());

    if( SOC[i] <= 0.0 )
    {
        SeB = SMALL;
    }
    else
    {
        SeB = (R*pow(I, 2) - 1.0*sign*I*T*dEdT)/V;
    }


{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rho*cp, T)
          - fvm::laplacian(Kappa, T)
          ==
            SeB
        );

        TEqn().relax();
        TEqn().solve();
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}
