    Info << "    Type: Natural Convection" << endl;

    const dimensionedVector& g = gFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    volScalarField& rhok = rhokFluid[i];
    rhok = 1.0 - beta*(T - TRef);

    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turb.divDevReff()
    );

    UEqn.relax();

    pZones.addResistance(UEqn);

    eqnResidual = solve
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
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);

