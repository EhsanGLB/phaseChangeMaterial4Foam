
if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/UEqnFC.H"

    for (int corr = 0; corr < nCorr; corr++)
    {
        #include "forcedConvection/pEqnFC.H"
    }
}


if(typeFluid[i] == "naturalConvection")
{
    if(PCMProperties[i].headerOk())
    {
        #include "PCM/UEqnPCM.H"

        for (int corr = 0; corr < nCorr; corr++)
        {
            #include "PCM/pEqnPCM.H"
        }
    }

    else
    {
        #include "naturalConvection/UEqnNC.H"

        for (int corr = 0; corr < nCorr; corr++)
        {
            #include "naturalConvection/pEqnNC.H"
        }
    }
}

