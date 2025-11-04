    const fvMesh& mesh = fluidRegions[i];

    volScalarField& T = TFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];


    volScalarField& rho = rhoFluid[i];
    volScalarField& Cp = CpFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volScalarField& mu = muFluid[i];
    volScalarField& beta = betaFluid[i];
    word& type = typeFluid[i];

    mesh.schemesDict().setFluxRequired(p.name());

   //-base fluid
    dimensionedScalar& rhof = rhofFluid[i];
    dimensionedScalar& Cpf = CpfFluid[i];
    dimensionedScalar& kappaf = kappafFluid[i];
    dimensionedScalar& muf = mufFluid[i];


    if (EPCMProperties[i].headerOk())
    {
        //-EPCM
        volScalarField& alphaEPCM = alphaEPCMFluid[i];
        volScalarField& rhoEPCM = rhoEPCMFluid[i];
        volScalarField& fEPCM = fEPCMFluid[i];
        volScalarField& moltenEPCM = moltenEPCMFluid[i];
        volScalarField& CpEPCM = CpEPCMFluid[i];
        word& alphaEPCMEq = alphaEPCMEqFluid[i];
        dimensionedScalar& omegaEPCM = omegaEPCMFluid[i];
        dimensionedScalar& lambdaEPCM = lambdaEPCMFluid[i];
        dimensionedScalar& kappaEPCM = kappaEPCMFluid[i];
        //-core
        dimensionedScalar& rhoc = rhocFluid[i];
        dimensionedScalar& Cpc = CpcFluid[i];
        dimensionedScalar& kappac = kappacFluid[i];
        dimensionedScalar& Tmc = TmcFluid[i];
        dimensionedScalar& Tmrc = TmrcFluid[i];
        dimensionedScalar& hmc = hmcFluid[i];
        dimensionedScalar& dc = dcFluid[i];
        volScalarField& Cpce = CpceFluid[i];
        //-shell
        dimensionedScalar& rhos = rhosFluid[i];
        dimensionedScalar& Cps = CpsFluid[i];
        dimensionedScalar& kappas = kappasFluid[i];
        dimensionedScalar& ds = dsFluid[i];
        //-EPCM
        dc = ds*pow( rhos/(rhos+omegaEPCM*rhoc), scalar(0.3333) );
        kappaEPCM = ds/( (dc/kappac) + ((ds-dc)/kappas) );

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

        rhoEPCM = ((1+omegaEPCM)*rhoc*rhos)/(rhos+omegaEPCM*rhoc);
        Cpce = Cpc+( (3.14/2)*(hmc/Tmrc-Cpc)*sin( (3.14/Tmrc)*(T-(Tmc-Tmrc/2)) ) )*fEPCM;
        CpEPCM = ((Cpce+omegaEPCM*Cps)*rhoc*rhos)/((rhos+omegaEPCM*rhoc)*rhoEPCM);

        rhoEPCM.correctBoundaryConditions();
        Cpce.correctBoundaryConditions();
        CpEPCM.correctBoundaryConditions();

        //-mixture
        rho = alphaEPCM*rhoEPCM+(1-alphaEPCM)*rhof;
        Cp = (alphaEPCM*rhoEPCM*CpEPCM+(1-alphaEPCM)*rhof*Cpf)/rho;
        kappa = ((kappaEPCM+(lambdaEPCM-1)*kappaf-(lambdaEPCM-1)*(kappaf-kappaEPCM)*alphaEPCM)/(kappaEPCM+(lambdaEPCM-1)*kappaf+(kappaf-kappaEPCM)*alphaEPCM))*kappaf;
        mu = (muf/pow((1.0-alphaEPCM),2.5));

        if(type == "naturalConvection")
        {
            dimensionedScalar& betaf = betafFluid[i];
            volScalarField& betaEPCM = betaEPCMFluid[i];
            dimensionedScalar& betac = betacFluid[i];
            dimensionedScalar& betas = betasFluid[i];
            betaEPCM = betac + ( (betas-betac)/2 )*(1 - omegaEPCM*rhos/rhoc);
            beta = (alphaEPCM*rhoEPCM*betaEPCM+(1-alphaEPCM)*rhof*betaf)/rho;
            beta.correctBoundaryConditions();
        }
    }
    else
    {
        rho = rhof;
        Cp = Cpf;
        kappa = kappaf;
        mu = muf;

        if(type == "naturalConvection")
        {
            dimensionedScalar& betaf = betafFluid[i];
            beta = betaf;
            beta.correctBoundaryConditions();
        }
    }

    rho.correctBoundaryConditions();
    Cp.correctBoundaryConditions();
    kappa.correctBoundaryConditions();
    mu.correctBoundaryConditions();



    incompressible::turbulenceModel& turb = turbulence[i];
    dimensionedScalar& Prt = PrtFluid[i];
    turb.correct();
    volScalarField& kappat = kappatFluid[i];
    kappat = (turb.nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    porousZones& pZones = pZonesFluid[i];
    volScalarField& porosity = porosityFluid[i];
    volScalarField& rhoPM = rhoPMFluid[i];
    volScalarField& CpPM = CpPMFluid[i];
    volScalarField& kappaPM = kappaPMFluid[i];
    volScalarField& hsPM = hsPMFluid[i];


    volScalarField& rhoEff = rhoEffFluid[i];
    volScalarField& CpEff = CpEffFluid[i];
    volScalarField& kappaEff = kappaEffFluid[i];


    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;
    rhoEff.correctBoundaryConditions();
    CpEff.correctBoundaryConditions();
    kappaEff.correctBoundaryConditions();


    dimensionedScalar& rh = rhFluid[i];

    if (resHeat[i].headerOk())
    {
        word& typerh = typerhFluid[i];
        dimensionedScalar& Rrh = RrhFluid[i];
        dimensionedScalar& Vmaxrh = VmaxrhFluid[i];
        dimensionedScalar& Imaxrh = ImaxrhFluid[i];
        dimensionedScalar& Tminrh = TminrhFluid[i];
        dimensionedScalar& Tmaxrh = TmaxrhFluid[i];
        dimensionedScalar Volrh = sum(fluidRegions[i].V());
        dimensionedScalar Pmaxrh = Vmaxrh*Imaxrh;
        dimensionedScalar PPrh ("PPrh", dimensionSet(1, 2, -3, 0, 0, 0, 0), SMALL);

        if (typerh == "none")
        {}

        else if (typerh == "power")
        {
            dimensionedScalar& Prh = PrhFluid[i];
            PPrh = Prh;
        }

        else if (typerh == "voltage")
        {
            dimensionedScalar& Vrh = VrhFluid[i];
            PPrh = Vrh*Vrh/Rrh;
        }

        else if (typerh == "current")
        {
            dimensionedScalar& Irh = IrhFluid[i];
            PPrh = Rrh*Irh*Irh;
        }

        else
        {
            WarningInFunction
                << "Resistive heating type '" << typerh << "' for region " << i
                << " is not recognized. Options are: none, power, voltage, current."
                << " Setting rh to 0." << endl;
        }

        if (Tminrh.value() > gMax(T) || gMax(T) > Tmaxrh.value())
        {
            rh = dimensionedScalar ("rh", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL);
            Info << "Resistive heating is inactive: " << Tminrh.value() << " < " << gMax(T) << " < " << Tmaxrh.value() << endl;
        }

        else
        {
            Info << "T max: " << gMax(T) << endl;
            rh = PPrh/Volrh;
            Info << "Resistive heating is active: " << PPrh.value() << " W or " << rh.value() << " W/m3" << endl;
        }
        
        if (PPrh.value() > Pmaxrh.value())
        {
            Info << "Warning: The generated power (" << PPrh.value() << ") is biger than the maximum power (" << Pmaxrh.value() << ")" << endl;
        }
    }



    radiation::radiationModel& rad = radiation[i];
