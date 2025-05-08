//  Pressure-velocity SIMPLE corrector

    p.storePrevIter();


if(typeFluid[i] == "forcedConvection")
{
    #include "forcedConvection/UEqnFC.H"
    #include "forcedConvection/pEqnFC.H"
    #include "forcedConvection/TEqnFC.H"

    if (alphaEPCMEqFluid[i] == "active")
    {
        #include "forcedConvection/alphaEqnFC.H"
    }
}


if(typeFluid[i] == "naturalConvection")
{
    if(PCMProperties[i].headerOk())
    {
        #include "PCM/UEqnPCM.H"
        #include "PCM/pEqnPCM.H"
        #include "PCM/TEqnPCM.H"

        if (alphaEPCMEqFluid[i] == "active")
        {
            #include "PCM/alphaEqnPCM.H"
        }
    }

    else
    {
        #include "naturalConvection/UEqnNC.H"
        #include "naturalConvection/pEqnNC.H"
        #include "naturalConvection/TEqnNC.H"

        if (alphaEPCMEqFluid[i] == "active")
        {
            #include "naturalConvection/alphaEqnNC.H"
        }
    }
}
