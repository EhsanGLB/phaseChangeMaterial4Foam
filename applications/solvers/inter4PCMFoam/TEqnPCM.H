//for (int oCorr=0; oCorr<6; oCorr++)
{
    //beta = (alpha1*rho1*beta1 + (scalar(1) - alpha1)*rho2*beta2)/rho;
    //volScalarField SePCM("SePCM", alpha1*SePCM1 + (scalar(1) - alpha1)*SePCM2);



    rhoCp = alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2;
    kappa = alpha1*kappa1 + (scalar(1) - alpha1)*kappa2;
    kappaEff = kappa;
    volScalarField SePCM1("SePCM1", hm1*4.0*exp(-pow(4.0*(T-Tm1)/Tmr1,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr1);
    volScalarField SePCM2("SePCM2", hm2*4.0*exp(-pow(4.0*(T-Tm2)/Tmr2,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr2);
    volScalarField& SePCMptr = SePCM();
    SePCMptr = rho1*SePCM1*alpha1 + rho2*SePCM2*(scalar(1) - alpha1);
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaEff, T)
      + SePCMptr*fvm::ddt(T)
      + SePCMptr*(U & fvc::grad(T))
      ==
        radiation->Ru()
      - fvm::Sp(4.0*radiation->Rp()*pow(T_,3), T)
      - radiation->Rp()*pow(T_,3)*(T_ - 4.0*T)
    );

    TEqn.relax();
    TEqn.solve();

    radiation->correct();
}
