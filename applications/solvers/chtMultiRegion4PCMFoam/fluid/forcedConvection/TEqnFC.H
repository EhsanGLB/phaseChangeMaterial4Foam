{
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
	fvm::ddt(rhoEff*CpEff, T)
      + fvm::div(phi*fvc::interpolate(rho*Cp), T)
      - fvm::laplacian(kappaEff, T)
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
