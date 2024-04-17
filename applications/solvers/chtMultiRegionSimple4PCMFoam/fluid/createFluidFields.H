    // Initialise fluid field pointer lists
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to pFluid\n" << endl;
        pFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
            )
        );
    }

//------------------- PCM4Foam --------------------//
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> CpFluid(fluidRegions.size());
    PtrList<volScalarField> kappaFluid(fluidRegions.size());
    PtrList<volScalarField> muFluid(fluidRegions.size());
    PtrList<volScalarField> betaFluid(fluidRegions.size());
    PtrList<singlePhaseTransportModel> laminarTransport(fluidRegions.size());
    PtrList<word> typeFluid(fluidRegions.size());

    scalarList initialMassFluid(fluidRegions.size());
    labelList pRefCellFluid(fluidRegions.size(),0);
    scalarList pRefValueFluid(fluidRegions.size(),0.0);

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rho", dimDensity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to CpFluid\n" << endl;
        CpFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Cp",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to kappaFluid\n" << endl;
        kappaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappa", dimThermalConductivity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to muFluid\n" << endl;
        muFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "mu",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to betaFluid\n" << endl;
        betaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "beta",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("beta", dimless/dimTemperature, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to laminarTransport\n" << endl;
        laminarTransport.set
        (
            i,
            new singlePhaseTransportModel
            (
                UFluid[i],
                phiFluid[i]
            )
        );

        typeFluid.set
        (
            i,
            new word (laminarTransport[i].lookup("type"))
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        setRefCell
        (
            pFluid[i],
            fluidRegions[i].solutionDict().subDict("SIMPLE"),
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

//------------------- PCM4Foam -> Base fluid --------------------//
    PtrList<dimensionedScalar> rhofFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappafFluid(fluidRegions.size());
    PtrList<dimensionedScalar> mufFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        rhofFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("rho"))
        );

        CpfFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Cp"))
        );

        kappafFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("kappa"))
        );

        mufFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("mu"))
        );
    }



    // forced convection
    PtrList<simpleControl> simpleFluid(fluidRegions.size());

    // natural convection and PCM
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betafFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TRefFluid(fluidRegions.size());
    PtrList<volScalarField> rhokFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
      // forced convection
      if(typeFluid[i] == "forcedConvection")
      {
        Info<< "    Adding to simpleFluid\n" << endl;
        simpleFluid.set
        (
            i,
            new simpleControl
            (
                fluidRegions[i]
            )
        );
      }


      // natural convection and PCM
      if(typeFluid[i] == "naturalConvection")
      {
        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        if(laminarTransport[i].found("beta"))
        {
            betafFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("beta"))
            );
        }

        if(laminarTransport[i].found("TRef"))
        {
            TRefFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("TRef"))
            );
        }

        Info<< "    Adding to rhokFluid\n" << endl;
        rhokFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhok", dimless, SMALL),
                "zeroGradient"
            )
        );
      }
    }


    // PCM
    PtrList<IOdictionary> PCMProperties(fluidRegions.size());
    PtrList<volScalarField> gammaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> hmfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmrfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CufFluid(fluidRegions.size());
    PtrList<dimensionedScalar> bfFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        PCMProperties.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "PCMProperties",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

      // PCM
      if(PCMProperties[i].headerOk())
      {
        Info<< "    Adding to gammaFluid\n" << endl;
        gammaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "gamma",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("gamma", dimless, SMALL),
                "zeroGradient"
            )
        );

        if(PCMProperties[i].found("hm"))
        {
            hmfFluid.set
            (
                i,
                new dimensionedScalar (PCMProperties[i].lookup("hm"))
            );
        }

        if(PCMProperties[i].found("Tm"))
        {
            TmfFluid.set
            (
                i,
                new dimensionedScalar (PCMProperties[i].lookup("Tm"))
            );
        }

        if(PCMProperties[i].found("Tmr"))
        {
            TmrfFluid.set
            (
                i,
                new dimensionedScalar (PCMProperties[i].lookup("Tmr"))
            );
        }

        if(PCMProperties[i].found("Cu"))
        {
            CufFluid.set
            (
                i,
                new dimensionedScalar (PCMProperties[i].lookup("Cu"))
            );
        }

        if(PCMProperties[i].found("b"))
        {
            bfFluid.set
            (
                i,
                new dimensionedScalar (PCMProperties[i].lookup("b"))
            );
        }
      }
    }

//------------------- PCM4Foam -> EPCM ------------------//
    PtrList<IOdictionary> EPCMProperties(fluidRegions.size());
    PtrList<volScalarField> alphaEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> rhoEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> fEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> moltenEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> CpEPCMFluid(fluidRegions.size());
    PtrList<volScalarField> betaEPCMFluid(fluidRegions.size());
    PtrList<word> alphaEPCMEqFluid(fluidRegions.size());
    PtrList<dimensionedScalar> omegaEPCMFluid(fluidRegions.size());
    PtrList<dimensionedScalar> lambdaEPCMFluid(fluidRegions.size());
    PtrList<dimensionedScalar> DBEPCMFluid(fluidRegions.size());
    PtrList<dimensionedScalar> DTEPCMFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappaEPCMFluid(fluidRegions.size());
    //-core
    PtrList<dictionary> coreSubDict(fluidRegions.size());
    PtrList<dimensionedScalar> rhocFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappacFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmrcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> hmcFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betacFluid(fluidRegions.size());
    PtrList<dimensionedScalar> dcFluid(fluidRegions.size());
    PtrList<volScalarField> CpceFluid(fluidRegions.size());
    //-shell
    PtrList<dictionary> shellSubDict(fluidRegions.size());
    PtrList<dimensionedScalar> rhosFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpsFluid(fluidRegions.size());
    PtrList<dimensionedScalar> kappasFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betasFluid(fluidRegions.size());
    PtrList<dimensionedScalar> dsFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
        EPCMProperties.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "EPCMProperties",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        //-EPCM
        if (EPCMProperties[i].headerOk())
        {
            coreSubDict.set
            (
                i,
                new dictionary
                (
                    EPCMProperties[i].subDict("core")
                )
            );

            shellSubDict.set
            (
                i,
                new dictionary
                (
                    EPCMProperties[i].subDict("shell")
                )
            );

            Info<< "    Adding to alphaEPCMFluid\n" << endl;
            alphaEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i]
                )
            );

            Info<< "    Adding to rhoEPCMFluid\n" << endl;
            rhoEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rhoEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("rhoEPCM", dimDensity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to fEPCMFluid\n" << endl;
            fEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "fEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("fEPCM", dimless, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to moltenEPCMFluid\n" << endl;
            moltenEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "moltenEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("moltenEPCM", dimless, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to CpEPCMFluid\n" << endl;
            CpEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "CpEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("CpEPCM", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to betaEPCMFluid\n" << endl;
            betaEPCMFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betaEPCM",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("betaEPCM", dimless/dimTemperature, SMALL),
                    "zeroGradient"
                )
            );

            alphaEPCMEqFluid.set
            (
                i,
                new word (EPCMProperties[i].lookup("alphaEPCMEq"))
            );

            omegaEPCMFluid.set
            (
                i,
                new dimensionedScalar (EPCMProperties[i].lookup("omega"))
            );

            lambdaEPCMFluid.set
            (
                i,
                new dimensionedScalar (EPCMProperties[i].lookup("lambda"))
            );

            if (alphaEPCMEqFluid[i] == "active")
            {
                DBEPCMFluid.set
                (
                    i,
                    new dimensionedScalar (EPCMProperties[i].lookup("DB"))
                );

                DTEPCMFluid.set
                (
                    i,
                    new dimensionedScalar (EPCMProperties[i].lookup("DT"))
                );
            }
            //-core
            rhocFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("rho"))
            );

            CpcFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("Cp"))
            );

            kappacFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("kappa"))
            );

            TmcFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("Tm"))
            );

            TmrcFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("Tmr"))
            );

            hmcFluid.set
            (
                i,
                new dimensionedScalar (coreSubDict[i].lookup("hm"))
            );

            if(coreSubDict[i].found("beta"))
            {
                betacFluid.set
                (
                    i,
                    new dimensionedScalar (coreSubDict[i].lookup("beta"))
                );
            }

            Info<< "    Adding to CpceFluid\n" << endl;
            CpceFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Cpce",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("Cpce", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            //-shell
            rhosFluid.set
            (
                i,
                new dimensionedScalar (shellSubDict[i].lookup("rho"))
            );

            CpsFluid.set
            (
                i,
                new dimensionedScalar (shellSubDict[i].lookup("Cp"))
            );

            kappasFluid.set
            (
                i,
                new dimensionedScalar (shellSubDict[i].lookup("kappa"))
            );

            if(shellSubDict[i].found("beta"))
            {
                betasFluid.set
                (
                    i,
                    new dimensionedScalar (shellSubDict[i].lookup("beta"))
                );
            }

            dsFluid.set
            (
                i,
                new dimensionedScalar (shellSubDict[i].lookup("d"))
            );

            kappaEPCMFluid.set
            (
                i,
                new dimensionedScalar ("kappaEPCM", dimThermalConductivity, SMALL)
            );

            dcFluid.set
            (
                i,
                new dimensionedScalar ("dc", dimLength, SMALL)
            );
        }
    }

//------------------- PCM4Foam -> Turbulence ------------------//
    PtrList<volScalarField> kappatFluid(fluidRegions.size());
    PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to kappat\n" << endl;
        kappatFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappat",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappat", dimThermalConductivity, SMALL)
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
                incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    laminarTransport[i]
                )
            ).ptr()
        );

        PrtFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Prt"))
        );
    }

//------------------- PCM4Foam -> Porous Media ------------------//
    PtrList<porousZones> pZonesFluid(fluidRegions.size());
    PtrList<volScalarField> porosityFluid(fluidRegions.size());
    PtrList<volScalarField> rhoPMFluid(fluidRegions.size());
    PtrList<volScalarField> CpPMFluid(fluidRegions.size());
    PtrList<volScalarField> kappaPMFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to pZones\n" << endl;
        pZonesFluid.set
        (
            i,
            new porousZones
            (
                fluidRegions[i]
            )
        );

        Info<< "    Adding to porosity\n" << endl;
        porosityFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "porosity",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("porosity", dimless, 1.0),
                "zeroGradient"
            )
        );

        Info<< "    Adding to rhoPM\n" << endl;
        rhoPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhoPM", dimDensity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to CpPM\n" << endl;
        CpPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "CpPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("CpPM", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        Info<< "    Adding to kappaPM\n" << endl;
        kappaPMFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaPM",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappaPM", dimThermalConductivity, SMALL),
                "zeroGradient"
            )
        );

        forAll( pZonesFluid[i], zoneI )
        {
            const label & zoneId( pZonesFluid[i][zoneI].zoneId() );
            const labelList & cells(fluidRegions[i].cellZones()[zoneId]);
            const scalar & zonePorosity( pZonesFluid[i][zoneI].porosity() );

            const dictionary& thermalPropertiesPM = pZonesFluid[i][zoneI].dict().subDict("thermalProperties");
            dimensionedScalar zoneRho(thermalPropertiesPM.lookup("rho"));
            dimensionedScalar zoneCp(thermalPropertiesPM.lookup("Cp"));
            dimensionedScalar zoneKappa(thermalPropertiesPM.lookup("kappa"));

            forAll( cells, cellI )
            {
                porosityFluid[i][cells[cellI]] = zonePorosity;
                rhoPMFluid[i][cells[cellI]] = zoneRho.value();
                CpPMFluid[i][cells[cellI]] = zoneCp.value();
                kappaPMFluid[i][cells[cellI]] = zoneKappa.value();
            }
        }
    }

//------------------- PCM4Foam -> Effective thermal properties ------------------//
    PtrList<volScalarField> kappaEffFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to kappaEff\n" << endl;
        kappaEffFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaEff",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("kappaEff", dimThermalConductivity, SMALL)
            )
        );
    }

//------------------- PCM4Foam -> Radiation ------------------//
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());


    forAll(fluidRegions, i)
    {
        Info<< "    Adding to radiation\n" << endl;
        radiation.set
        (
            i,
            radiation::radiationModel::New
            (
                TFluid[i]
            )
        );

    }
