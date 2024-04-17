{
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
	fvm::ddt(rhoEff*CpEff, T)
      + rho*Cp*fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
      ==
        rad.Ru()
      - fvm::Sp(4.0*rad.Rp()*pow(T_,3), T)
      - rad.Rp()*pow(T_,3)*(T_ - 4.0*T)
    );

    TEqn.relax();
    TEqn.solve();

    rad.correct();
}
