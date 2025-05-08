    Info << "    Type: Forced Convection" << endl;
    simpleControl& simplec = simpleFluid[i];

    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turb.divDevReff()
    );

     pZones.addResistance(UEqn());

    // Momentum solution
    eqnResidual = solve(relax(UEqn()) == -fvc::grad(p)).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
