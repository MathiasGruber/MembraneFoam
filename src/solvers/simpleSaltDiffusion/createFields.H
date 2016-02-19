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