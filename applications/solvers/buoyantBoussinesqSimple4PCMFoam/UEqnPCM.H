    // Solve the momentum equation

    volScalarField& GAMMA = gamma();
    GAMMA = 0.5*erf(4.0*(T-Tmf)/Tmrf)+scalar(0.5);
    volScalarField SmPCM("SmPCM", Cu*pow((1.0-GAMMA), 2.0)/(pow(GAMMA, 3.0)+b));

    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
      + (1/rho)*fvm::SuSp(SmPCM, U)
    );

    UEqn.relax();

    pZones.addResistance(UEqn);

    solve
    (
        UEqn
        ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(rhok)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
            )
        )
    );


