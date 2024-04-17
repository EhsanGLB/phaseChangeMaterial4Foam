    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
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

#   include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rho*gh
    );


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, pimple.dict(), pdRefCell, pdRefValue);
    mesh.schemesDict().setFluxRequired(pd.name());

    scalar pRefValue = 0.0;

    if (pd.needReference())
    {
        pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pdRefCell)
        );
    }

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );


//-------------------PCM4Foam-->energy------------------//
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


    const word& phase1Name = twoPhaseProperties.phase1Name();
    const dictionary& phase1SubDict = twoPhaseProperties.subDict(phase1Name);
    dimensionedScalar Cp1(phase1SubDict.lookup("Cp"));
    dimensionedScalar kappa1(phase1SubDict.lookup("kappa"));
    dimensionedScalar mu1(phase1SubDict.lookup("mu"));
    dimensionedScalar beta1(phase1SubDict.lookup("beta"));
    dimensionedScalar Prt1(phase1SubDict.lookup("Prt"));

    const word& phase2Name = twoPhaseProperties.phase2Name();
    const dictionary& phase2SubDict = twoPhaseProperties.subDict(phase2Name);
    dimensionedScalar Cp2(phase1SubDict.lookup("Cp"));
    dimensionedScalar kappa2(phase1SubDict.lookup("kappa"));
    dimensionedScalar mu2(phase1SubDict.lookup("mu"));
    dimensionedScalar beta2(phase2SubDict.lookup("beta"));
    dimensionedScalar Prt2(phase2SubDict.lookup("Prt"));

    dimensionedScalar TRef(twoPhaseProperties.lookup("TRef"));

    Info<< "Reading field Cp\n" << endl;
    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Reading field kappa\n" << endl;
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Reading field kappaEff\n" << endl;
    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaEff", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Reading field beta\n" << endl;
    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("beta", dimless/dimTemperature, SMALL),
        "zeroGradient"
    );


    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2,
        alpha1.boundaryField().types()
    );
    rhoCp.oldTime();

    surfaceScalarField rhoCpPhi
    (
        IOobject
        (
            "rhoCpPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoPhi*Cp1
    );


    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - (alpha1*beta1 + (scalar(1) - alpha1)*beta2)*(T - TRef)
    );

    volScalarField rhoBeta
    (
        IOobject
        (
            "rhoBeta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho1*beta1*alpha1 + rho2*beta2*(scalar(1) - alpha1)
    );


//-------------------PCM4Foam-->PCM------------------//
    IOdictionary PCMProperties
    (
        IOobject
        (
            "PCMProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    autoPtr<volScalarField> gamma1(nullptr);
    autoPtr<volScalarField> SmPCM(nullptr);
    autoPtr<volScalarField> SePCM(nullptr);
    dimensionedScalar hm1("hm1", dimEnergy/dimMass, SMALL);
    dimensionedScalar Tm1("Tm1", dimTemperature, SMALL);
    dimensionedScalar Tmr1("Tmr1", dimTemperature, SMALL);
    dimensionedScalar Cu1("Cu1", dimMass/dimVolume/dimTime, SMALL);
    dimensionedScalar b1("b1", dimless, SMALL);


    if (PCMProperties.headerOk())
    {
        Info<< "Reading field gamma1\n" << endl;
        gamma1.set
        (
            new volScalarField
            (
                IOobject
                (
                    "gamma1",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("gamma1", dimless, SMALL),
                "zeroGradient"
            )
        );

        Info<< "Reading field SmPCM\n" << endl;
        SmPCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "SmPCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("SmPCM", dimensionSet(1,-3,-1,0,0,0,0), SMALL),
                "zeroGradient"
            )
        );

        Info<< "Reading field SePCM\n" << endl;
        SePCM.set
        (
            new volScalarField
            (
                IOobject
                (
                    "SePCM",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("SePCM", dimDensity*dimEnergy/dimMass/dimTemperature, SMALL),
                "zeroGradient"
            )
        );

        const dictionary& phase1SubDictPCM = PCMProperties.subDict(phase1Name);
        hm1 = dimensionedScalar(phase1SubDictPCM.lookup("hm"));
        Tm1 = dimensionedScalar(phase1SubDictPCM.lookup("Tm"));
        Tmr1 = dimensionedScalar(phase1SubDictPCM.lookup("Tmr"));
        Cu1 = dimensionedScalar(phase1SubDictPCM.lookup("Cu"));
        b1 = dimensionedScalar(phase1SubDictPCM.lookup("b"));
    }


    autoPtr<volScalarField> gamma2(nullptr);
    dimensionedScalar hm2("hm2", dimEnergy/dimMass, SMALL);
    dimensionedScalar Tm2("Tm2", dimTemperature, SMALL);
    dimensionedScalar Tmr2("Tmr2", dimTemperature, SMALL);
    dimensionedScalar Cu2("Cu2", dimMass/dimVolume/dimTime, SMALL);
    dimensionedScalar b2("b2", dimless, SMALL);


    if (PCMProperties.headerOk())
    {
        Info<< "Reading field gamma2\n" << endl;
        gamma2.set
        (
            new volScalarField
            (
                IOobject
                (
                    "gamma2",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("gamma2", dimless, SMALL),
                "zeroGradient"
            )
        );

        const dictionary& phase2SubDictPCM = PCMProperties.subDict(phase2Name);
        hm2 = dimensionedScalar(phase2SubDictPCM.lookup("hm"));
        Tm2 = dimensionedScalar(phase2SubDictPCM.lookup("Tm"));
        Tmr2 = dimensionedScalar(phase2SubDictPCM.lookup("Tmr"));
        Cu2 = dimensionedScalar(phase2SubDictPCM.lookup("Cu"));
        b2 = dimensionedScalar(phase2SubDictPCM.lookup("b"));
    }

//-------------------PCM4Foam-->EPCM------------------//




//-------------------PCM4Foam-->Radiation------------------//
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );
