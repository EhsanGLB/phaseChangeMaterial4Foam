// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

// Convection-diffusion matrix
fvVectorMatrix HUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevReff()
);

pZones.addResistance(HUEqn);

if (pimple.momentumPredictor())
{
    // Solve momentum predictor
    solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
}
else
{
    // Explicit update
    U = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p))/(HUEqn.A() + ddtUEqn.A());
    U.correctBoundaryConditions();
}
