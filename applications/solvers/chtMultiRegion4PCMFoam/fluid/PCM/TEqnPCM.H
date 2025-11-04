{
    dimensionedScalar& hmf = hmfFluid[i];
    dimensionedScalar& Tmf = TmfFluid[i];
    dimensionedScalar& Tmrf = TmrfFluid[i];

    volScalarField SePCM("SePCM", hmf*4.0*exp(-pow(4.0*(T-Tmf)/Tmrf,2))/Foam::sqrt(mathematicalConstant::pi)/Tmrf);
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
	fvm::ddt(rhoEff*CpEff, T)
      + fvm::div(phi*fvc::interpolate(rho*Cp), T)
      - fvm::laplacian(kappaEff, T)
      + rho*porosity*SePCM*fvm::ddt(T)
      + rho*porosity*SePCM*(U & fvc::grad(T))
      ==
        rad.Ru()
      - fvm::Sp(4.0*rad.Rp()*pow(T_,3), T)
      - rad.Rp()*pow(T_,3)*(T_ - 4.0*T)
      + hsPM
      + rh
    );

    TEqn.relax();
    TEqn.solve();

    rad.correct();
}
