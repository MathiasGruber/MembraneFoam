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

    dimensionedScalar pi_mACoeff(transportProperties.lookup("pi_mACoeff"));
    dimensionedScalar mu0(transportProperties.lookup("mu0"));
    dimensionedScalar mu_mACoeff(transportProperties.lookup("mu_mACoeff"));
    dimensionedScalar D_AB_Min(transportProperties.lookup("D_AB_Min"));
    dimensionedScalar D_AB_Coeff(transportProperties.lookup("D_AB_Coeff"));
    dimensionedScalar D_AB_mACoeff(transportProperties.lookup("D_AB_mACoeff"));
    dimensionedScalar rho0(transportProperties.lookup("rho0"));
    dimensionedScalar rho_mACoeff(transportProperties.lookup("rho_mACoeff"));

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


    // surfaceScalarField mW
    // (
    //     IOobject
    //     (
    //         "Mw",
    //         runTime.timeName(),
    //         mesh,
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