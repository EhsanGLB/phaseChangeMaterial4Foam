{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();

    surfaceScalarField phiU
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, U, phi)
    );

    phi = phiU + rUAf*fvc::interpolate(rhok)*(g & mesh.Sf());

    //while (pimple.correctNonOrthogonal())
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p) == fvc::div(phi)
        );

        //pEqn.solve(mesh.solutionDict().solver(p.select(pimple.finalInnerIter())));

        if
        (
            oCorr == nOuterCorr-1
         && corr == nCorr - 1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
        }
        else
        {
            pEqn.solve(mesh.solutionDict().solver(p.name()));
        }


        //if (pimple.finalNonOrthogonalIter())
        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    #include "incompressibleContinuityErrs.H"
}
