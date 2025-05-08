{
    rhoCp = alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2;
    kappa = alpha1*kappa1 + (scalar(1) - alpha1)*kappa2;
    kappaEff = kappa;
    volScalarField T_ = T;


    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaEff, T)
      ==
        radiation->Ru()
      - fvm::Sp(4.0*radiation->Rp()*pow(T_,3), T)
      - radiation->Rp()*pow(T_,3)*(T_ - 4.0*T)
    );

    TEqn.relax();
    TEqn.solve();

    radiation->correct();
}
