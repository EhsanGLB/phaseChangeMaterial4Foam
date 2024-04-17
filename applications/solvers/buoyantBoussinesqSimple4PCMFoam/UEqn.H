    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
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
