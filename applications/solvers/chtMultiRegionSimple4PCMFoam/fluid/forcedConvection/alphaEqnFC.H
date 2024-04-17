{
    volScalarField& alphaEPCM = alphaEPCMFluid[i];
    dimensionedScalar& DBEPCM = DBEPCMFluid[i];
    dimensionedScalar& DTEPCM = DTEPCMFluid[i];

    fvScalarMatrix alphaEqn
    (
        fvm::div(phi, alphaEPCM)
      - fvm::laplacian(DBEPCM, alphaEPCM)
      - fvc::laplacian((DTEPCM/T), T)
    );

    alphaEqn.relax();
    alphaEqn.solve();
}

