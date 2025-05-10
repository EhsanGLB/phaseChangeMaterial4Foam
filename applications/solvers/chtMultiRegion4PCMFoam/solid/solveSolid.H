if (battery[i].headerOk())
{
    #include "battery.H"
}

else if (thermoelectricGenerator[i].headerOk())
{
    #include "thermoelectricGenerator.H"
}

else if(metalHydride[i].headerOk())
{
    #include "metalHydride.H"
}

else if (EPCMPropertiesSolid[i].headerOk())
{
    #include "EPCM.H"
}

else
{
    #include "TEqn.H"
}

