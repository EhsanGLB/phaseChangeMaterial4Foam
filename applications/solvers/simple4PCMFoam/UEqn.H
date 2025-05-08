    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    pZones.addResistance(UEqn());

    // Momentum solution
    solve(relax(UEqn()) == -fvc::grad(p));
