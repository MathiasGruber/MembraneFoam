/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pisoSaltTransport
    based on twoLiquidMixingFoam and rhoPisoFoam

Description
    Solver for soluble salt transport in a fluid. 

    The solver assumes the process is isothermal and uses a weakly 
    compressible solution technique. The density is solely a function of the 
    solute mass fraction and is calculated explicitly.

    The current model does not include a turbulence model as the standard
    compressible models require more sophisticated energy treatments. The 
    incompressible turbulence models are not used as the density variation is
    significant enough that a weakly-compressible formulation is required.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readControls.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "m_AInitialContinuity.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Outer-corrector loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            #include "UEqn.H"

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "m_AEqn.H"
                #include "pEqn.H"
            }
        }

        // Calculate density
        rho = rho0 * (1.0 + rho_mACoeff*m_A);

        runTime.write();

        #include "m_AContinuity.H"

        // Time since last timestep
        scalar timeSinceLast = runTime.clockTimeIncrement();
        Info<< "Time since last iteration = " << timeSinceLast << " s" << endl;

        // Show the execution speed
        Info<< "Simulation Speed = "
            <<  runTime.deltaTValue() / ( timeSinceLast / 86400 )
            << " s / day" << endl;

        // Show the execution time
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
