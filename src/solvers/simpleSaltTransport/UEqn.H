    volScalarField mu
    (
        "mu",
        mu0 * (1 + mu_mACoeff * m_A)
    );

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(mu, U)
      - fvc::div(mu*dev2(fvc::grad(U)().T()))
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    // if (simple.momentumPredictor())
    // {
    //     solve(UEqn == -fvc::grad(p));

    //     fvOptions.correct(U);
    // }