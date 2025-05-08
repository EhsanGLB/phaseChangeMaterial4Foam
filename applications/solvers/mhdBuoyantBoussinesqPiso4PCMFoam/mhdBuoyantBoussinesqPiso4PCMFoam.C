/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantBoussinesqPisoFoam

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pisoControl.H"
#include "radiationModel.H"
#include "porousZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
#       include "readBPISOControls.H"

        Info << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        while (piso.loop())
        {
#           include "updateProperties.H"

            if (PCMProperties.headerOk())
            {
#               include "UEqnPCM.H"

                while (piso.correct())
                {
#                   include "pEqn.H"
                }

#               include "TEqnPCM.H"

                if (alphaEPCMEq == "active")
                {
#                   include "alphaEqn.H"
                }

		// --- B-PISO loop
#		include "pBEqn.H"
            }

            else
            {
#               include "UEqn.H"

                while (piso.correct())
                {
#                   include "pEqn.H"
                }

#               include "TEqn.H"

                if (alphaEPCMEq == "active")
                {
#                   include "alphaEqn.H"
                }

		// --- B-PISO loop
#		include "pBEqn.H"
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
