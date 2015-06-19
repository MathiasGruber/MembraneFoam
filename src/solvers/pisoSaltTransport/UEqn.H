    volScalarField mu
    (
        "mu",
        mu0 + mu_mACoeff*m_A
    );
 
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      - fvm::laplacian(mu, U)
      - fvc::div(mu*dev2(fvc::grad(U)().T()))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn 
        == 
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }
