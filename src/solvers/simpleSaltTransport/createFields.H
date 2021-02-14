Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field m_A\n" << endl;
    volScalarField m_A
    (
        IOobject
        (
            "m_A",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );
    dimensionedScalar pi_mACoeff("pi_mACoeff", dimless, transportProperties);
    dimensionedScalar mu0("mu0", dimDynamicViscosity, transportProperties);
    dimensionedScalar mu_mACoeff("mu_mACoeff", dimless, transportProperties);
    dimensionedScalar D_AB_Min("D_AB_Min", dimViscosity, transportProperties);
    dimensionedScalar D_AB_Coeff("D_AB_Coeff", dimViscosity, transportProperties);
    dimensionedScalar D_AB_mACoeff("D_AB_mACoeff", dimless, transportProperties);
    dimensionedScalar rho0("rho0", dimDensity, transportProperties);
    dimensionedScalar rho_mACoeff("rho_mACoeff", dimless, transportProperties);

    Info << "Creating rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho0 * (1.0 + rho_mACoeff*m_A)
    );
    rho.oldTime();

    // Info << "Creating Jw\n" << endl;
    // // surfaceScalarField Jw
    // volVectorField Jw
    // (
    //     IOobject
    //     (
    //         "Jw",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     // mesh
    // );
    // Jw.oldTime();

    // surfaceScalarField Jw
    // (
    //     IOobject
    //     (
    //         "Jw",
    //         runTime.timeName(),
    //         mesh,
    //         // patch().boundaryMesh().mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::interpolate(U) & mesh.Sf()
    // );


    #include "compressibleCreatePhi.H"

    Info << "Flow is assumed to be laminar\nTurbulence is currently unsupported in this solver\n" << endl;

    Info<< "Calculating field (g.h)f\n" << endl;
    surfaceScalarField ghf = surfaceScalarField("ghf", g & mesh.Cf());

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
