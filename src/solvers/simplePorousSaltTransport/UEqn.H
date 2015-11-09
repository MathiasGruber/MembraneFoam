    // Update the viscosity field
    mu = mu0 + mu_mACoeff*m_A;

    // Construct the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
          fvm::div(phi, U)
          - fvm::laplacian(mu, U)
          - fvc::div(mu*dev2(fvc::grad(U)().T()))
    );
    UEqn().relax();

    // Update the mu_laplacian field
    mu_laplacianU = fvc::laplacian( mu, U );

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal
    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn().A();
        pZones.addResistance(UEqn(), tTU());
        trTU = inv(tTU());
        trTU().rename("rAU");

        fvOptions.constrain(UEqn());

        volVectorField gradp(fvc::grad(p));

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn().H() - gradp);
        }
        U.correctBoundaryConditions();

        fvOptions.correct(U);
    }
    else
    {
        pZones.addResistance(UEqn());
        fvOptions.constrain(UEqn());
        solve(UEqn() == -fvc::grad(p));
        fvOptions.correct(U);

        trAU = 1.0/UEqn().A();
        trAU().rename("rAU");
    }
