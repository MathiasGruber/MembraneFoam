    // Load transport properties dictionary
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

    // Load relevant transport properties
    dimensionedScalar pi_mACoeff(transportProperties.lookup("pi_mACoeff"));
    dimensionedScalar mu0(transportProperties.lookup("mu0"));
    dimensionedScalar mu_mACoeff(transportProperties.lookup("mu_mACoeff"));
    dimensionedScalar D_AB_Min(transportProperties.lookup("D_AB_Min"));
    dimensionedScalar D_AB_Coeff(transportProperties.lookup("D_AB_Coeff"));
    dimensionedScalar D_AB_mACoeff(transportProperties.lookup("D_AB_mACoeff"));
    dimensionedScalar rho0(transportProperties.lookup("rho0"));
    dimensionedScalar rho_mACoeff(transportProperties.lookup("rho_mACoeff"));

    // Load / Create fields
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

    Info<< "Reading field Diff_ratio\n" << endl;
    volScalarField Diff_ratio
    (
        IOobject
        (
            "Diff_ratio",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field inVesicle\n" << endl;
    volScalarField inVesicle
    (
        IOobject
        (
            "inVesicle",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    Info << "Creating mu\n" << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mu0 + mu_mACoeff*m_A
    );

    Info << "Creating mu*laplacianU\n" << endl;
    volVectorField mu_laplacianU
    (
        IOobject
        (
            "mu_laplacianU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mu * fvc::laplacian( U )
    );

    // Create Identity matrix
    volVectorField UnityVectorField
    (
        IOobject
        (
            "UnityVectorField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("UnityVectorField",dimensionSet(0,0,0,0,0,0,0),Foam::vector(1,1,1) )
    );

    Info << "Creating rP\n" << endl;
    volScalarField rP
    (
        IOobject
        (
            "rP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
   

    #include "compressibleCreatePhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);
