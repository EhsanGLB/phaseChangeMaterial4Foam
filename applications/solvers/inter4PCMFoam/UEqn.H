
    rhoBeta = rho1*beta1*alpha1 + rho2*beta2*(scalar(1) - alpha1);

    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
      ==
      - rhoBeta*(T-TRef)*g
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
               (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                )*mesh.magSf()
            )

          /*+ fvc::reconstruct
            (
               (
                    fvc::interpolate(rhoBeta)*(g & mesh.Sf())
                )
            )*/
        );
    }

