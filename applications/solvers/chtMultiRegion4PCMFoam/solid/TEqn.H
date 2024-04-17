{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rho*cp, T)
          - fvm::laplacian(Kappa, T)
          ==
            hs
          + SeB
          + SeTEG
        );

        TEqn().relax();
        TEqn().solve();

        if (battery[i].headerOk())
        {
            #include "battery.H"
        }

        if (thermoelectricGenerator[i].headerOk())
        {
            #include "thermoelectricGenerator.H"
        }
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}
