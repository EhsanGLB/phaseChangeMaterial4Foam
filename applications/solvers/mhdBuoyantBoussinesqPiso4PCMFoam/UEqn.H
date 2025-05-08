    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
      - fvc::div(phiB, 2.0*DBU*B)
      + fvc::grad(DBU*magSqr(B))
    );

    UEqn.relax();

    pZones.addResistance(UEqn);

    if (piso.momentumPredictor())
    {
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
    }
