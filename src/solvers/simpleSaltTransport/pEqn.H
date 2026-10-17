{
    // Calculate new density
    rho = rho0 * (1.0 + rho_mACoeff * m_A);
    
    Info<< "rho max/min : "
        << max(rho).value() << " "
        << min(rho).value() << endl;
    
    Info<< "pressure max/min : "
        << max(p).value() << " "
        << min(p).value() << endl;
    
    Info<< "m_A max/min : "
        << max(m_A).value() << " "
        << min(m_A).value() << endl;
    
    Info<< "U max/min : "
        << max(U).value() << " "
        << min(U).value() << endl;

    // Calculate rAU
    volScalarField rAU(1.0/UEqn().A());
    
    // Interpolate rAU to faces
    tmp<surfaceScalarField> temp1 = fvc::interpolate(rAU);
    surfaceScalarField rAUf = temp1();
    
    // Calculate phiU
    // Preserve prescribed membrane/inlet velocities on the pressure predictor.
    volVectorField HbyA(constrainHbyA(rAU*UEqn().H(), U, p));
    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(rho) *
        (
            fvc::interpolate(HbyA) & mesh.Sf()
        )
    );
    
    // Calculate Phi
    phi = phiU - fvc::interpolate(rho)*ghf*fvc::snGrad(rho)*rAUf*mesh.magSf();

    // Solve the equation
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rAU, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }
    

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Update velocity BCs
    // Use relaxed pressure in the velocity corrector. Reconstructing from
    // the unrelaxed pressure flux made the old pressure relaxation ineffective.
    U = HbyA - rAU*
    (
        fvc::grad(p)
      + fvc::reconstruct(ghf*fvc::snGrad(rho)*mesh.magSf())
    );
    U.correctBoundaryConditions();

}
