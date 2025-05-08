Info << "    Type: Forced Convection" << endl;
volScalarField& rAU = rAUFluid[i];
pimpleControl& pimplec = pimpleFluid[i];


// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

// Convection-diffusion matrix
fvVectorMatrix HUEqn
(
    fvm::div(phi, U)
  + turb.divDevReff()
);

pZones.addResistance(HUEqn);

if (momentumPredictor)
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
