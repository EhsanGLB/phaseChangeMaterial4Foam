
{
    p.boundaryField().updateCoeffs();

    // Prepare clean 1/a_p without time derivative and under-relaxation
    // contribution
    rAU = 1.0/HUEqn.A();

    // Calculate U from convection-diffusion matrix
    U = rAU*HUEqn.H();

    // Consistently calculate flux
    pimplec.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

    // Global flux balance
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    //while (pimple.correctNonOrthogonal())
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAU)/pimplec.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
        );

        //pEqn.setReference(pRefCell, pRefValue);
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

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    // Explicitly relax pressure for momentum corrector except for last
    // corrector
    /*if (!pimple.finalIter())
    {
        p.relax();
    }*/

#   include "incompressibleContinuityErrs.H"

    // Consistently reconstruct velocity after pressure equation. Note: flux is
    // made relative inside the function
    pimplec.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
}
