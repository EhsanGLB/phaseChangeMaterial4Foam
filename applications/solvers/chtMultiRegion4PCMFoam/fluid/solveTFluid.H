
if(typeFluid[i] == "forcedConvection")
{
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
        #include "PCM/TEqnPCM.H"

        if (alphaEPCMEqFluid[i] == "active")
        {
            #include "PCM/alphaEqnPCM.H"
        }
    }

    else
    {
        #include "naturalConvection/TEqnNC.H"

        if (alphaEPCMEqFluid[i] == "active")
        {
            #include "naturalConvection/alphaEqnNC.H"
        }
    }
}


