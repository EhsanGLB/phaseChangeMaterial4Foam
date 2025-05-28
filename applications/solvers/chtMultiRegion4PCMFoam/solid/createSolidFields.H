    // Initialise solid field pointer lists
    PtrList<volScalarField> rhos(solidRegions.size());
    PtrList<volScalarField> cps(solidRegions.size());
    PtrList<volScalarField> rhosCps(solidRegions.size());
    PtrList<volScalarField> Kappas(solidRegions.size());
    PtrList<volScalarField> Ts(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to rhos\n" << endl;
        rhos.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to cps\n" << endl;
        cps.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "cp",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        rhosCps.set
        (
            i,
            new volScalarField("rhosCps", rhos[i]*cps[i])
        );

        Info<< "    Adding to Kappas\n" << endl;
        Kappas.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Kappa",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to Ts\n" << endl;
        Ts.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
    }

//-------------------PCM4Foam-->Battery------------------//
    PtrList<IOdictionary> battery(solidRegions.size());
    PtrList<word> stateBattery(solidRegions.size());
    PtrList<scalar> CBattery(solidRegions.size());
    PtrList<scalar> ISOCBattery(solidRegions.size());
    PtrList<dimensionedScalar> RBattery(solidRegions.size());
    PtrList<dimensionedScalar> IBattery(solidRegions.size());
    PtrList<dimensionedScalar> dEdTBattery(solidRegions.size());
    scalarList SOC(solidRegions.size());
    PtrList<volScalarField> SeBattery(solidRegions.size());


    forAll(solidRegions, i)
    {
        battery.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "battery",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        SeBattery.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "SeB",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar("SeB", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL),
                "zeroGradient"
            )
        );

        if(battery[i].found("state"))
        {
            stateBattery.set
            (
                i,
                new word (battery[i].lookup("state"))
            );
        }

        if(battery[i].found("C"))
        {
            CBattery.set
            (
                i,
                new scalar (readScalar(battery[i].lookup("C")))
            );
        }

        if(battery[i].found("ISOC"))
        {
            ISOCBattery.set
            (
                i,
                new scalar (readScalar(battery[i].lookup("ISOC")))
            );

            SOC[i] = ISOCBattery[i];
        }

        if(battery[i].found("R"))
        {
            RBattery.set
            (
                i,
                new dimensionedScalar (battery[i].lookup("R"))
            );
        }

        if(battery[i].found("I"))
        {
            IBattery.set
            (
                i,
                new dimensionedScalar (battery[i].lookup("I"))
            );
        }

        if(battery[i].found("dEdT"))
        {
            dEdTBattery.set
            (
                i,
                new dimensionedScalar (battery[i].lookup("dEdT"))
            );
        }
    }

    scalar tOld = 0.0;

//-------------------PCM4Foam-->Heat source------------------//
    PtrList<IOdictionary> heatSource(solidRegions.size());
    PtrList<dimensionedScalar> hsSolid(solidRegions.size());

    forAll(solidRegions, i)
    {
        heatSource.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "heatSource",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        if(heatSource[i].found("hs"))
        {
            hsSolid.set
            (
                i,
                new dimensionedScalar (heatSource[i].lookup("hs"))
            );
        }
        else
        {
            hsSolid.set
            (
                i,
                new dimensionedScalar ("hs", dimThermalConductivity*dimTemperature/dimArea, SMALL)
            );
        }
    }

//-------------------PCM4Foam-->Thermoelectric generator------------------//
    PtrList<IOdictionary> thermoelectricGenerator(solidRegions.size());
    PtrList<dimensionedScalar> ZTTEG(solidRegions.size());
    PtrList<dimensionedScalar> SeTEGG(solidRegions.size());

    forAll(solidRegions, i)
    {
        thermoelectricGenerator.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "thermoelectricGenerator",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        SeTEGG.set
        (
            i,
            new dimensionedScalar ("SeTEG", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL)
        );

        if(thermoelectricGenerator[i].found("ZT"))
        {
            ZTTEG.set
            (
                i,
                new dimensionedScalar (thermoelectricGenerator[i].lookup("ZT"))
            );
        }
    }

//-------------------PCM4Foam-->EPCM------------------//
    PtrList<IOdictionary> EPCMPropertiesSolid(solidRegions.size());
    PtrList<volScalarField> alphaEPCMSolid(solidRegions.size());
    PtrList<volScalarField> rhoEPCMSolid(solidRegions.size());
    PtrList<volScalarField> fEPCMSolid(solidRegions.size());
    PtrList<volScalarField> moltenEPCMSolid(solidRegions.size());
    PtrList<volScalarField> CpEPCMSolid(solidRegions.size());
    PtrList<volScalarField> rhomSolid(solidRegions.size());
    PtrList<volScalarField> CpmSolid(solidRegions.size());
    PtrList<volScalarField> kappamSolid(solidRegions.size());

    PtrList<dimensionedScalar> omegaEPCMSolid(solidRegions.size());
    PtrList<dimensionedScalar> lambdaEPCMSolid(solidRegions.size());
    PtrList<dimensionedScalar> kappaEPCMSolid(solidRegions.size());
    //-core
    PtrList<dictionary> coreSubDictSolid(solidRegions.size());
    PtrList<dimensionedScalar> rhocSolid(solidRegions.size());
    PtrList<dimensionedScalar> CpcSolid(solidRegions.size());
    PtrList<dimensionedScalar> kappacSolid(solidRegions.size());
    PtrList<dimensionedScalar> TmcSolid(solidRegions.size());
    PtrList<dimensionedScalar> TmrcSolid(solidRegions.size());
    PtrList<dimensionedScalar> hmcSolid(solidRegions.size());
    PtrList<dimensionedScalar> dcSolid(solidRegions.size());
    PtrList<volScalarField> CpceSolid(solidRegions.size());
    //-shell
    PtrList<dictionary> shellSubDictSolid(solidRegions.size());
    PtrList<dimensionedScalar> rhoshSolid(solidRegions.size());
    PtrList<dimensionedScalar> CpshSolid(solidRegions.size());
    PtrList<dimensionedScalar> kappashSolid(solidRegions.size());
    PtrList<dimensionedScalar> dshSolid(solidRegions.size());


    forAll(solidRegions, i)
    {
        EPCMPropertiesSolid.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "EPCMProperties",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        //-EPCM
        if (EPCMPropertiesSolid[i].headerOk())
        {
            coreSubDictSolid.set
            (
                i,
                new dictionary
                (
                    EPCMPropertiesSolid[i].subDict("core")
                )
            );

            shellSubDictSolid.set
            (
                i,
                new dictionary
                (
                    EPCMPropertiesSolid[i].subDict("shell")
                )
            );

            Info<< "    Adding to alphaEPCMSolid\n" << endl;
            alphaEPCMSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "alphaEPCMs",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i]
                )
            );

            Info<< "    Adding to rhoEPCMSolid\n" << endl;
            rhoEPCMSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rhoEPCMs",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("rhoEPCMs", dimDensity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to fEPCMSolid\n" << endl;
            fEPCMSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "fEPCMs",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("fEPCMs", dimless, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to moltenEPCMSolid\n" << endl;
            moltenEPCMSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "moltenEPCMs",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("moltenEPCMs", dimless, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to CpEPCMSolid\n" << endl;
            CpEPCMSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "CpEPCMs",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("CpEPCMs", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to rhomSolid\n" << endl;
            rhomSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rhoms",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("rhoms", dimDensity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to CpmSolid\n" << endl;
            CpmSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Cpms",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("Cpms", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to kappamSolid\n" << endl;
            kappamSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "kappams",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("kappams", dimThermalConductivity, SMALL),
                    "zeroGradient"
                )
            );

            omegaEPCMSolid.set
            (
                i,
                new dimensionedScalar (EPCMPropertiesSolid[i].lookup("omega"))
            );

            lambdaEPCMSolid.set
            (
                i,
                new dimensionedScalar (EPCMPropertiesSolid[i].lookup("lambda"))
            );

            //-core
            rhocSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("rho"))
            );

            CpcSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("Cp"))
            );

            kappacSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("kappa"))
            );

            TmcSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("Tm"))
            );

            TmrcSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("Tmr"))
            );

            hmcSolid.set
            (
                i,
                new dimensionedScalar (coreSubDictSolid[i].lookup("hm"))
            );

            Info<< "    Adding to CpceSolid\n" << endl;
            CpceSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Cpces",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("Cpces", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            //-shell
            rhoshSolid.set
            (
                i,
                new dimensionedScalar (shellSubDictSolid[i].lookup("rho"))
            );

            CpshSolid.set
            (
                i,
                new dimensionedScalar (shellSubDictSolid[i].lookup("Cp"))
            );

            kappashSolid.set
            (
                i,
                new dimensionedScalar (shellSubDictSolid[i].lookup("kappa"))
            );

            dshSolid.set
            (
                i,
                new dimensionedScalar (shellSubDictSolid[i].lookup("d"))
            );

            kappaEPCMSolid.set
            (
                i,
                new dimensionedScalar ("kappaEPCM", dimThermalConductivity, SMALL)
            );

            dcSolid.set
            (
                i,
                new dimensionedScalar ("dc", dimLength, SMALL)
            );
        }
    }

//-------------------PCM4Foam-->Metal hydride storage------------------//
    PtrList<IOdictionary> metalHydride(solidRegions.size());
    PtrList<dictionary> gasSubDict(solidRegions.size());
    PtrList<dictionary> chargingPressureSubDict(solidRegions.size());
    PtrList<volScalarField> rhoEffMetalHydride(solidRegions.size());
    PtrList<volScalarField> CpEffMetalHydride(solidRegions.size());
    PtrList<volScalarField> KappaEffMetalHydride(solidRegions.size());
    PtrList<volScalarField> pEqMetalHydride(solidRegions.size());
    PtrList<volScalarField> SmMHMetalHydride(solidRegions.size());
    PtrList<volScalarField> SeMHMetalHydride(solidRegions.size());


    PtrList<dimensionedScalar> deltaHMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> deltaSMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> CaMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> EaMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> rhosMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> pRefMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> RMetalHydride(solidRegions.size());


    PtrList<dimensionedScalar> rhogMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> CpgMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> kappagMetalHydride(solidRegions.size());
    PtrList<dimensionedScalar> MgMetalHydride(solidRegions.size());


    PtrList<scalar> aMetalHydride(solidRegions.size());
    PtrList<scalar> bMetalHydride(solidRegions.size());

    forAll(solidRegions, i)
    {
        metalHydride.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "metalHydride",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        if (metalHydride[i].headerOk())
        {

            gasSubDict.set
            (
                i,
                new dictionary
                (
                    metalHydride[i].subDict("gas")
                )
            );

            chargingPressureSubDict.set
            (
                i,
                new dictionary
                (
                    metalHydride[i].subDict("chargingPressure")
                )
            );

            Info<< "    Adding to rhoEff\n" << endl;
            rhoEffMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rhoEff",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("rhoEff", dimDensity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to CpEff\n" << endl;
            CpEffMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "CpEff",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("CpEff", dimSpecificHeatCapacity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to KappaEff\n" << endl;
            KappaEffMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "KappaEff",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("KappaEff", dimThermalConductivity, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to pEq\n" << endl;
            pEqMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "pEq",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("pEq", dimPressure, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to SmMH\n" << endl;
            SmMHMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "SmMH",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("SmMH", dimDensity/dimTime, SMALL),
                    "zeroGradient"
                )
            );

            Info<< "    Adding to SmMH\n" << endl;
            SeMHMetalHydride.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "SeMH",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("SeMH", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL),
                    "zeroGradient"
                )
            );

            deltaHMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("deltaH"))
            );

            deltaSMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("deltaS"))
            );

            CaMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("Ca"))
            );

            EaMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("Ea"))
            );

            rhosMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("rhos"))
            );

            pRefMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("pRef"))
            );

            RMetalHydride.set
            (
                i,
                new dimensionedScalar (metalHydride[i].lookup("R"))
            );

            rhogMetalHydride.set
            (
                i,
                new dimensionedScalar (gasSubDict[i].lookup("rho"))
            );

            CpgMetalHydride.set
            (
                i,
                new dimensionedScalar (gasSubDict[i].lookup("Cp"))
            );

            kappagMetalHydride.set
            (
                i,
                new dimensionedScalar (gasSubDict[i].lookup("kappa"))
            );

            MgMetalHydride.set
            (
                i,
                new dimensionedScalar (gasSubDict[i].lookup("M"))
            );

            aMetalHydride.set
            (
                i,
                new scalar (readScalar(chargingPressureSubDict[i].lookup("a")))
            );

            bMetalHydride.set
            (
                i,
                new scalar (readScalar(chargingPressureSubDict[i].lookup("b")))
            );


        }
    }

//-------------------PCM4Foam-->Porous Media------------------//
    PtrList<porousZones> pZonesSolid(solidRegions.size());
    PtrList<volScalarField> porositySolid(solidRegions.size());

    forAll(solidRegions, i)
    {
        Info<< "    Adding to pZones\n" << endl;
        pZonesSolid.set
        (
            i,
            new porousZones
            (
                solidRegions[i]
            )
        );

        porositySolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "porosity",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar("porosity", dimless, 1.0),
                "zeroGradient"
            )
        );

        forAll( pZonesSolid[i], zoneI )
        {
            const label & zoneId( pZonesSolid[i][zoneI].zoneId() );
            const labelList & cells(solidRegions[i].cellZones()[zoneId]);
            const scalar & zonePorosity( pZonesSolid[i][zoneI].porosity() );

            forAll( cells, cellI )
            {
                porositySolid[i][cells[cellI]] = zonePorosity;
            }
        }

        //porositySolid[i].write();
    }

