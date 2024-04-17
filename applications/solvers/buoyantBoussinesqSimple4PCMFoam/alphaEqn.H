{
    volScalarField& alphaEPCMptr = alphaEPCM();

    fvScalarMatrix alphaEqn
    (
        fvm::div(phi, alphaEPCMptr)
      + fvm::SuSp(-fvc::div(phi), alphaEPCMptr)
      - fvm::laplacian(DBEPCM, alphaEPCMptr)
      - fvc::laplacian((DTEPCM/T), T)
    );

    alphaEqn.relax();
    alphaEqn.solve();
}

