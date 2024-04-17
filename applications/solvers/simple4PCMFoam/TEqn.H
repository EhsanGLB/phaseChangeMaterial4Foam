{
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
        rho*Cp*fvm::div(phi, T)
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
