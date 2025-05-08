fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff()
  - fvc::div(phiB, 2.0*DBU*B)
  + fvc::grad(DBU*magSqr(B))
);

pZones.addResistance(UEqn);
solve(UEqn == -fvc::grad(p));
