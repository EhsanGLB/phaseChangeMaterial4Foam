    Info << "    Type: PCM" << endl;

    volScalarField& GAMMA1 = gamma1();
    volScalarField& GAMMA2 = gamma2();
    GAMMA1 = 0.5*erf(4.0*(T-Tm1)/Tmr1)+scalar(0.5);
    GAMMA2 = 0.5*erf(4.0*(T-Tm2)/Tmr2)+scalar(0.5);
    rhoBeta = rho1*beta1*alpha1 + rho2*beta2*(scalar(1) - alpha1);

    volScalarField SmPCM1("SmPCM1", Cu1*pow((1.0-GAMMA1), 2.0)/(pow(GAMMA1, 3.0)+b1));
    volScalarField SmPCM2("SmPCM2", Cu2*pow((1.0-GAMMA2), 2.0)/(pow(GAMMA2, 3.0)+b2));
    volScalarField& SmPCMptr = SmPCM();
    SmPCMptr = alpha1*SmPCM1 + (scalar(1) - alpha1)*SmPCM2;

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
      + fvm::SuSp(SmPCMptr, U)
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
        );
    }

