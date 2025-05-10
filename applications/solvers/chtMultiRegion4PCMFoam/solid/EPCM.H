//-EPCM
volScalarField& alphaEPCM = alphaEPCMSolid[i];
volScalarField& rhoEPCM = rhoEPCMSolid[i];
volScalarField& fEPCM = fEPCMSolid[i];
volScalarField& moltenEPCM = moltenEPCMSolid[i];
volScalarField& CpEPCM = CpEPCMSolid[i];
volScalarField& rhoms = rhomSolid[i];
volScalarField& Cpms = CpmSolid[i];
volScalarField& kappams = kappamSolid[i];
dimensionedScalar& omegaEPCM = omegaEPCMSolid[i];
dimensionedScalar& lambdaEPCM = lambdaEPCMSolid[i];
dimensionedScalar& kappaEPCM = kappaEPCMSolid[i];
//-core
dimensionedScalar& rhoc = rhocSolid[i];
dimensionedScalar& Cpc = CpcSolid[i];
dimensionedScalar& kappac = kappacSolid[i];
dimensionedScalar& Tmc = TmcSolid[i];
dimensionedScalar& Tmrc = TmrcSolid[i];
dimensionedScalar& hmc = hmcSolid[i];
dimensionedScalar& dc = dcSolid[i];
volScalarField& Cpce = CpceSolid[i];
//-shell
dimensionedScalar& rhosh = rhoshSolid[i];
dimensionedScalar& Cpsh = CpshSolid[i];
dimensionedScalar& kappash = kappashSolid[i];
dimensionedScalar& dsh = dshSolid[i];
//-EPCM
dc = dsh*pow( rhosh/(rhosh+omegaEPCM*rhoc), scalar(0.3333) );
kappaEPCM = dsh/( (dc/kappac) + ((dsh-dc)/kappash) );

forAll(fEPCM, j)
{
    if (T[j] < scalar(Tmc.value()-Tmrc.value()/2) )
    {
        fEPCM[j] = scalar(0.0);
        moltenEPCM[j] = scalar(0.0);
    }
    else if ( scalar(Tmc.value()-Tmrc.value()/2) <= T[j] && T[j] <= scalar(Tmc.value()+Tmrc.value()/2) )
    {
        fEPCM[j] = scalar(1.0);
        moltenEPCM[j] = scalar(0.5);
    }
    else if (T[j] > scalar(Tmc.value()+Tmrc.value()/2) )
    {
        fEPCM[j] = scalar(0.0);
        moltenEPCM[j] = scalar(1.0);
    }
}

fEPCM.correctBoundaryConditions();
moltenEPCM.correctBoundaryConditions();

rhoEPCM = ((1+omegaEPCM)*rhoc*rhosh)/(rhosh+omegaEPCM*rhoc);
Cpce = Cpc+( (3.14/2)*(hmc/Tmrc-Cpc)*sin( (3.14/Tmrc)*(T-(Tmc-Tmrc/2)) ) )*fEPCM;
CpEPCM = ((Cpce+omegaEPCM*Cpsh)*rhoc*rhosh)/((rhosh+omegaEPCM*rhoc)*rhoEPCM);

rhoEPCM.correctBoundaryConditions();
Cpce.correctBoundaryConditions();
CpEPCM.correctBoundaryConditions();

//-mixture
#include "setRegionSolidFields.H"
rhoms = alphaEPCM*rhoEPCM+(1-alphaEPCM)*rho;
Cpms = (alphaEPCM*rhoEPCM*CpEPCM+(1-alphaEPCM)*rho*cp)/rhoms;
kappams = ((kappaEPCM+(lambdaEPCM-1)*Kappa-(lambdaEPCM-1)*(Kappa-kappaEPCM)*alphaEPCM)/(kappaEPCM+(lambdaEPCM-1)*Kappa+(Kappa-kappaEPCM)*alphaEPCM))*Kappa;


{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rhoms*Cpms, T)
          - fvm::laplacian(kappams, T)
        );

        TEqn().relax();
        TEqn().solve();
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}
