
    if (EPCMProperties.headerOk())
    {
        volScalarField& alphaEPCMptr = alphaEPCM();
        volScalarField& fEPCMptr = fEPCM();
        volScalarField& moltenEPCMptr = moltenEPCM();
        volScalarField& CpEPCMptr = CpEPCM();
        volScalarField& Cpceptr = Cpce();

        forAll(fEPCMptr, j)
        {
            if (T[j] < scalar(Tmc.value()-Tmrc.value()/2) )
            {
                fEPCMptr[j] = scalar(0.0);
                moltenEPCMptr[j] = scalar(0.0);
            }
            else if ( scalar(Tmc.value()-Tmrc.value()/2) <= T[j] && T[j] <= scalar(Tmc.value()+Tmrc.value()/2) )
            {
                fEPCMptr[j] = scalar(1.0);
                moltenEPCMptr[j] = scalar(0.5);
            }
            else if (T[j] > scalar(Tmc.value()+Tmrc.value()/2) )
            {
                fEPCMptr[j] = scalar(0.0);
                moltenEPCMptr[j] = scalar(1.0);
            }
        }

        fEPCMptr.correctBoundaryConditions();
        moltenEPCMptr.correctBoundaryConditions();


        Cpceptr = Cpc+( (3.14/2)*(hmc/Tmrc-Cpc)*sin( (3.14/Tmrc)*(T-(Tmc-Tmrc/2)) ) )*fEPCMptr;
        CpEPCMptr = ((Cpceptr+omegaEPCM*Cps)*rhoc*rhos)/((rhos+omegaEPCM*rhoc)*rhoEPCM);

        Cpceptr.correctBoundaryConditions();
        CpEPCMptr.correctBoundaryConditions();

        //-mixture
        rho = alphaEPCMptr*rhoEPCM+(1-alphaEPCMptr)*rhof;
        Cp = (alphaEPCMptr*rhoEPCM*CpEPCMptr+(1-alphaEPCMptr)*rhof*Cpf)/rho;
        kappa = ((kappaEPCM+(lambdaEPCM-1)*kappaf-(lambdaEPCM-1)*(kappaf-kappaEPCM)*alphaEPCMptr)/(kappaEPCM+(lambdaEPCM-1)*kappaf+(kappaf-kappaEPCM)*alphaEPCMptr))*kappaf;
        mu = (muf/pow((1.0-alphaEPCMptr),2.5));
    }
    else
    {
        rho = rhof;
        Cp = Cpf;
        kappa = kappaf;
        mu = muf;
    }

    rho.correctBoundaryConditions();
    Cp.correctBoundaryConditions();
    kappa.correctBoundaryConditions();
    mu.correctBoundaryConditions();

    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;
    rhoEff.correctBoundaryConditions();
    CpEff.correctBoundaryConditions();
    kappaEff.correctBoundaryConditions();


    rhoR = rhoRf;
    DB = rhoR/mu0;
    DBU = 1.0/(2.0*mu0*rho);
