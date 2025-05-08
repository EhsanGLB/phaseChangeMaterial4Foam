{
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoEff*CpEff, T)
      + fvm::div(phi*fvc::interpolate(rho*Cp), T)
      - fvm::laplacian(kappaEff, T)
      - jh*(1/rhoR)*pow(mag(U^B), 2)
      ==
        radiation->Ru()
      - fvm::Sp(4.0*radiation->Rp()*pow(T_,3), T)
      - radiation->Rp()*pow(T_,3)*(T_ - 4.0*T)
      + hsPM
    );

    TEqn.relax();
    TEqn.solve();

    radiation->correct();
}
