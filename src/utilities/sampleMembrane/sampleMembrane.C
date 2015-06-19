/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    sampleMembrane

Description
    Samples the p, U and m_A field of a membrane patch and prints to post processing file

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::validArgs.append("patchName");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const word patchName = args[1];

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Get the fields
        IOobject VelocityHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        IOobject saltHeader
        (
            "m_A",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        IOobject PressureHeader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        IOobject ShearHeader
        (
            "wallShearStress",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check fields exists
        if (VelocityHeader.headerOk() && saltHeader.headerOk() && PressureHeader.headerOk())
        {
            mesh.readUpdate();

            // Check that the patch is there
            const label patchI = mesh.boundaryMesh().findPatchID(patchName);
            if (patchI < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }

            // Give patch area
            Info<< "    Area vector of patch "
                << patchName << '[' << patchI << ']' << " = "
                << gSum(mesh.Sf().boundaryField()[patchI]) << endl;
            Info<< "    Area magnitude of patch "
                << patchName << '[' << patchI << ']' << " = "
                << gSum(mesh.magSf().boundaryField()[patchI]) << endl;

            // Write out face centers of patch
            IOField<vector> cfOut
            (
                IOobject
                (
                    patchName+"_coordinates",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh.Cf().boundaryField()[ patchI ]
            );
            cfOut.write();
            
            // Write out face area of patch
            IOField<scalar> sfOut
            (
                IOobject
                (
                    patchName+"_face_areas",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh.magSf().boundaryField()[ patchI ]
            );
            sfOut.write();
            
            // Write out the other fields
            Info<< "    Reading m_A" << endl;
            volScalarField saltField(saltHeader, mesh);
            
            Info<< "    Reading U" << endl;
            volVectorField velocityField(VelocityHeader, mesh);
            
            Info<< "    Reading p" << endl;
            volScalarField pressureField(PressureHeader, mesh);
            
            // Write out the values
            IOField<scalar> saltOut
            (
                IOobject
                (
                    saltField.name()+"_"+patchName+"_values",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                saltField.boundaryField()[ patchI ]
            );
            saltOut.write();
            
            IOField<scalar> pressureOut
            (
                IOobject
                (
                    pressureField.name()+"_"+patchName+"_values",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pressureField.boundaryField()[ patchI ]
            );
            pressureOut.write();
            
            IOField<vector> velocityOut
            (
                IOobject
                (
                    velocityField.name()+"_"+patchName+"_values",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                velocityField.boundaryField()[ patchI ]
            );
            velocityOut.write();
            
            // Also output shear stresses
            if( ShearHeader.headerOk() ){
                Info<< "    Reading wallShearStress" << endl;
                volVectorField shearField(ShearHeader, mesh);
                IOField<vector> shearOut
                (
                    IOobject
                    (
                        shearField.name()+"_"+patchName+"_values",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    shearField.boundaryField()[ patchI ]
                );
                shearOut.write();
            }

            
        }
        else
        {
            Info<< "    Could not find all fields. p, U and m_A required. " << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
