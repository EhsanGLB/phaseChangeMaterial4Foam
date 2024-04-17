    fvMesh& mesh = solidRegions[i];

    volScalarField& rho = rhos[i];
    volScalarField& cp = cps[i];
    volScalarField& Kappa = Kappas[i];
    volScalarField& T = Ts[i];

    dimensionedScalar& hs = hsSolid[i];

    volScalarField& SeB = SeBattery[i];
    dimensionedScalar& SeTEG = SeTEGG[i];

