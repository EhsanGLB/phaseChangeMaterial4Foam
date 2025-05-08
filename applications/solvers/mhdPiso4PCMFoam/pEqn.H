volScalarField rUA = 1.0/UEqn.A();

U = rUA*UEqn.H();

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rUA, U, phi);

while (piso.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();
