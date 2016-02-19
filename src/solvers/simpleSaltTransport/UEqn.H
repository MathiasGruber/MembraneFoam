    volScalarField mu
    (
        "mu",
        mu0 * (1 + mu_mACoeff * m_A)
    );

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(mu, U)
      - fvc::div(mu*dev2(fvc::grad(U)().T()))
    );

    UEqn().relax();