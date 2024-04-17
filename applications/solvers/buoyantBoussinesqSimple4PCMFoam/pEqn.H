{
    volScalarField rUA("rUA", 1.0/UEqn.A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, U, phi)
    );

    phi = phiU + rUAf*fvc::interpolate(rhok)*(g & mesh.Sf());

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p) == fvc::div(phi)
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
