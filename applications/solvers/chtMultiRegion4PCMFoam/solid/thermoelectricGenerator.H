
    dimensionedScalar& ZT = ZTTEG[i];
    wordList TEGPatcheNames (thermoelectricGenerator[i].lookup("TEGPatcheNames"));
    surfaceScalarField gradT = fvc::snGrad(T);

    label patch0ID = solidRegions[i].boundaryMesh().findPatchID(TEGPatcheNames[0]);
    scalar patch0Heat = sum(Kappa.boundaryField()[patch0ID]*gradT.boundaryField()[patch0ID]*mesh.magSf().boundaryField()[patch0ID]);
    scalar patch0TAvg = sum(T.boundaryField()[patch0ID]*mesh.magSf().boundaryField()[patch0ID])/sum(solidRegions[i].magSf().boundaryField()[patch0ID]);

    label patch1ID = solidRegions[i].boundaryMesh().findPatchID(TEGPatcheNames[1]);
    scalar patch1Heat = sum(Kappa.boundaryField()[patch1ID]*gradT.boundaryField()[patch1ID]*mesh.magSf().boundaryField()[patch1ID]);
    scalar patch1TAvg = sum(T.boundaryField()[patch1ID]*solidRegions[i].magSf().boundaryField()[patch1ID])/sum(solidRegions[i].magSf().boundaryField()[patch1ID]);

    scalar etaTEG = ( (Foam::sqrt(1.0+ZT.value()) - scalar(1.0)) / (Foam::sqrt(1.0+ZT.value()) + (min(patch0TAvg, patch1TAvg)/max(patch0TAvg, patch1TAvg)) ) ) * (mag(patch0TAvg - patch1TAvg)/max(patch0TAvg, patch1TAvg));

    scalar heatTransfer = mag(patch1Heat-patch0Heat);
    dimensionedScalar P("P", dimensionSet(1, -1, -3, 0, 0, 0, 0),  (etaTEG*heatTransfer)/sum(solidRegions[i].V()).value());
    SeTEG = -1*P;
