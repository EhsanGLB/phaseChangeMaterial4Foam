    p.boundaryField().updateCoeffs();

    volScalarField rAU = 1.0/UEqn().A();
    U = rAU*UEqn().H();
    UEqn.clear();

    // Calculate under-relaxation consistent flux
    simplec.calcSteadyConsistentFlux(phi, U);
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (simplec.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAU)/simplec.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        eqnResidual = pEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

        if (simplec.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Reconstruct consistent velocity after pressure equation
    simplec.reconstructSteadyVelocity(U, rAU, p);
