    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//------------------- PCM4Foam ------------------//
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Constructing field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field Cp\n" << endl;
    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field kappa\n" << endl;
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field mu\n" << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field rhoR\n" << endl;
    volScalarField rhoR
    (
        IOobject
        (
            "rhoR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoR", dimensionSet(1, 3, -3, 0, 0, -2, 0), SMALL),
        "zeroGradient"
    );

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    singlePhaseTransportModel laminarTransport(U, phi);

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());


    Info<< "Reading transport Properties\n" << endl;
    dimensionedScalar rhof(laminarTransport.lookup("rho"));
    dimensionedScalar Cpf(laminarTransport.lookup("Cp"));
    dimensionedScalar kappaf(laminarTransport.lookup("kappa"));
    dimensionedScalar muf(laminarTransport.lookup("mu"));
    dimensionedScalar rhoRf(laminarTransport.lookup("rhoR"));

//------------------- PCM4Foam -> EPCM ------------------//
    IOdictionary EPCMProperties
    (
        IOobject
        (
            "EPCMProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<volScalarField> alphaEPCM(nullptr);
    autoPtr<volScalarField> fEPCM(nullptr);
    autoPtr<volScalarField> moltenEPCM(nullptr);
    autoPtr<volScalarField> CpEPCM(nullptr);

    word alphaEPCMEq;
    dimensionedScalar omegaEPCM("omegaEPCM", dimless, SMALL);
    dimensionedScalar lambdaEPCM("lambdaEPCM", dimless, SMALL);
    dimensionedScalar DBEPCM("DBEPCM", dimArea/dimTime, SMALL);
    dimensionedScalar DTEPCM("DTEPCM", dimArea/dimTime, SMALL);
    dimensionedScalar rhoEPCM("rhoEPCM", dimDensity, SMALL);
    dimensionedScalar kappaEPCM("kappaEPCM", dimThermalConductivity, SMALL);

//- EPCM's core
    dimensionedScalar rhoc("rhoc", dimDensity, SMALL);
    dimensionedScalar Cpc("Cpc", dimSpecificHeatCapacity, SMALL);
    dimensionedScalar kappac("kappac", dimThermalConductivity, SMALL);
    dimensionedScalar hmc("hmc", dimEnergy/dimMass, SMALL);
    dimensionedScalar Tmc("Tmc", dimTemperature, SMALL);
    dimensionedScalar Tmrc("Tmrc", dimTemperature, SMALL);
    dimensionedScalar dc("dc", dimLength, SMALL);
    autoPtr<volScalarField> Cpce(nullptr);

//- EPCM's shell
    dimensionedScalar rhos("rhos", dimDensity, SMALL);
    dimensionedScalar Cps("Cps", dimSpecificHeatCapacity, SMALL);
    dimensionedScalar kappas("kappas", dimThermalConductivity, SMALL);
    dimensionedScalar ds("ds", dimLength, SMALL);


    if (EPCMProperties.headerOk())
    {
        //- EPCM
        Info<< "Constructing field alphaEPCM\n" << endl;
        alphaEPCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "alphaEPCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "Constructing field fEPCM\n" << endl;
        fEPCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "fEPCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("fEPCM", dimless, SMALL),
                "zeroGradient"
            )
        );

        Info<< "Constructing field moltenEPCM\n" << endl;
        moltenEPCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "moltenEPCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("moltenEPCM", dimless, SMALL),
                "zeroGradient"
            )
        );

        Info<< "Constructing field CpEPCM\n" << endl;
        CpEPCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "CpEPCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("CpEPCM", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );


        alphaEPCMEq = word(EPCMProperties.lookup("alphaEPCMEq"));
        omegaEPCM = dimensionedScalar(EPCMProperties.lookup("omega"));
        lambdaEPCM = dimensionedScalar(EPCMProperties.lookup("lambda"));
        if (alphaEPCMEq == "active")
        {
            DBEPCM = dimensionedScalar(EPCMProperties.lookup("DB"));
            DTEPCM = dimensionedScalar(EPCMProperties.lookup("DT"));
        }

        //- EPCM's core
        dictionary coreSubDict = EPCMProperties.subDict("core");
        rhoc = dimensionedScalar(coreSubDict.lookup("rho"));
        Cpc = dimensionedScalar(coreSubDict.lookup("Cp"));
        kappac = dimensionedScalar(coreSubDict.lookup("kappa"));
        hmc = dimensionedScalar(coreSubDict.lookup("hm"));
        Tmc = dimensionedScalar(coreSubDict.lookup("Tm"));
        Tmrc = dimensionedScalar(coreSubDict.lookup("Tmr"));

        Info<< "Constructing field Cpce\n" << endl;
        Cpce.set
        (
            new volScalarField
            (
                IOobject
                (
                    "Cpce",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("Cpce", dimSpecificHeatCapacity, SMALL),
                "zeroGradient"
            )
        );

        //- EPCM's shell
        dictionary shellDict = EPCMProperties.subDict("shell");
        rhos = dimensionedScalar(shellDict.lookup("rho"));
        Cps = dimensionedScalar(shellDict.lookup("Cp"));
        kappas = dimensionedScalar(shellDict.lookup("kappa"));
        ds = dimensionedScalar(shellDict.lookup("d"));

        //- EPCM
        dc = ds*pow(rhos/(rhos+omegaEPCM*rhoc), scalar(0.3333));
        rhoEPCM = ((1+omegaEPCM)*rhoc*rhos)/(rhos+omegaEPCM*rhoc);
        kappaEPCM = ds/((dc/kappac)+((ds-dc)/kappas));
    }

//------------------- PCM4Foam -> Turbulence ------------------//
    Info<< "Constructing field kappat\n" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappat", dimThermalConductivity, SMALL)
    );

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

//------------------- PCM4Foam -> Porous Media ------------------//
    porousZones pZones(mesh);

    IOdictionary porousZonesDict
    (
        IOobject
        (
            "porousZones",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("porosity", dimless, 1.0),
        "zeroGradient"
    );

    volScalarField rhoPM
    (
        IOobject
        (
            "rhoPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPM", dimDensity, SMALL),
        "zeroGradient"
    );

    volScalarField CpPM
    (
        IOobject
        (
            "CpPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CpPM", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    volScalarField kappaPM
    (
        IOobject
        (
            "kappaPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaPM", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    volScalarField hsPM
    (
        IOobject
        (
            "hsPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("hsPM", dimensionSet(1, -1, -3, 0, 0, 0, 0), SMALL),
        "zeroGradient"
    );

    forAll( pZones, zoneI )
    {
        const label & zoneId( pZones[zoneI].zoneId() );
        const labelList & cells(mesh.cellZones()[zoneId]);
        const scalar & zonePorosity( pZones[zoneI].porosity() );

        const dictionary& thermalPropertiesPM = pZones[zoneI].dict().subDict("thermalProperties");
        dimensionedScalar zoneRho(thermalPropertiesPM.lookup("rho"));
        dimensionedScalar zoneCp(thermalPropertiesPM.lookup("Cp"));
        dimensionedScalar zoneKappa(thermalPropertiesPM.lookup("kappa"));
        dimensionedScalar zoneHs(thermalPropertiesPM.lookup("hs"));

        forAll( cells, cellI )
        {
            porosity[cells[cellI]] = zonePorosity;
            rhoPM[cells[cellI]] = zoneRho.value();
            CpPM[cells[cellI]] = zoneCp.value();
            kappaPM[cells[cellI]] = zoneKappa.value();
            hsPM[cells[cellI]] = zoneHs.value();
        }
    }

    rhoPM.correctBoundaryConditions();
    CpPM.correctBoundaryConditions();
    kappaPM.correctBoundaryConditions();
    hsPM.correctBoundaryConditions();

//------------------- PCM4Foam -> Effective thermal properties ------------------//
    Info<< "Constructing field rhoEff\n" << endl;
    volScalarField rhoEff
    (
        IOobject
        (
            "rhoEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoEff", dimDensity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field CpEff\n" << endl;
    volScalarField CpEff
    (
        IOobject
        (
            "CpEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("CpEff", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field kappaEff\n" << endl;
    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaEff", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

//------------------- PCM4Foam -> Radiation ------------------//
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );

//------------------- PCM4Foam -> Magnetic Field ------------------//

    IOdictionary mhdProperties
    (
        IOobject
        (
            "mhdProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar mu0(mhdProperties.lookup("mu0"));
    word jouleHeating = word(mhdProperties.lookup("jouleHeating"));
    scalar jh = 0.0;
    if (jouleHeating == "active")
    {
        jh = 1.0;
    }

    Info<< "Reading field pB\n" << endl;
    volScalarField pB
    (
        IOobject
        (
            "pB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field B\n" << endl;
    volVectorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


#   include "createPhiB.H"

    volScalarField DB("DB", rhoR/mu0);
    volScalarField DBU("DBU", 1.0/(2.0*mu0*rho));

    mesh.schemesDict().setFluxRequired(pB.name());
