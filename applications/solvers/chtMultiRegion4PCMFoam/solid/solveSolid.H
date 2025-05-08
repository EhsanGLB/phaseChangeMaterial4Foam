if (battery[i].headerOk())
{
    #include "battery.H"
}

if (thermoelectricGenerator[i].headerOk())
{
    #include "thermoelectricGenerator.H"
}

if(metalHydride[i].headerOk())
{
    #include "metalHydride.H"
}

else
{
    #include "TEqn.H"
}

