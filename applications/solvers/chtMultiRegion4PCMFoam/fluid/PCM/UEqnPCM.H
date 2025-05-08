    Info << "    Type: PCM" << endl;

    const dimensionedVector& g = gFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    volScalarField& rhok = rhokFluid[i];
    volScalarField& gamma = gammaFluid[i];
    dimensionedScalar& hmf = hmfFluid[i];
    dimensionedScalar& Tmf = TmfFluid[i];
    dimensionedScalar& Tmrf = TmrfFluid[i];
    dimensionedScalar& Cuf = CufFluid[i];
    dimensionedScalar& bf = bfFluid[i];

    // Solve the Momentum equation
    rhok = 1.0 - beta*(T - TRef);
    gamma = 0.5 * erf( 4.0 * ( T - Tmf ) / Tmrf ) + scalar(0.5);
    volScalarField SmPCM("SmPCM", Cuf * pow( ( 1 - gamma ) , 2 ) / ( pow( gamma , 3 ) + bf));
    Info << "max(p): " << max(p) << endl;


    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turb.divDevReff()
      + (1/rho)*fvm::SuSp(SmPCM, U)
    );

    UEqn().relax();

    pZones.addResistance(UEqn());

    if (momentumPredictor)
    {
         solve
         (
            UEqn()
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rhok)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
