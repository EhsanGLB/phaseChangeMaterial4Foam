{
    volScalarField& alphaEPCMptr = alphaEPCM();

    fvScalarMatrix alphaEqn
    (
        fvm::ddt(alphaEPCMptr)
      + fvm::div(phi, alphaEPCMptr)
      - fvm::laplacian(DBEPCM, alphaEPCMptr)
      - fvc::laplacian((DTEPCM/T), T)
    );

    alphaEqn.relax();
    alphaEqn.solve();
}

