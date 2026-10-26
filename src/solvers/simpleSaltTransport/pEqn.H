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

    // Consistent pressure defect correction for compatible Laplacians with
    // linear coefficient interpolation. Include density and gravity in the
    // cell correction and use the discrete matrix flux on faces.
    tmp<volScalarField> consistentRAU(nullptr);
    if (simple.consistent())
    {
        forAll(p.boundaryField(), patchi)
        {
            const word kind = p.boundaryField()[patchi].type();
            if (kind != "fixedValue" && kind != "zeroGradient"
             && kind != "empty" && kind != "processor"
             && kind != "symmetryPlane" && kind != "symmetry")
                FatalErrorInFunction << "Consistent correction requires a linear "
                    << "fixedValue/zeroGradient/empty/processor/symmetry pressure patch; "
                    << "found " << kind << exit(FatalError);
        }
        const volScalarField denominator(1.0/rAU - UEqn().H1());
        const scalar minDen = gMin(denominator.primitiveField());
        const scalar maxDen = gMax(denominator.primitiveField());
        if (!(minDen > 0) || !std::isfinite(minDen) || !std::isfinite(maxDen))
            FatalErrorInFunction << "Non-positive/non-finite consistent momentum "
                << "diagonal " << minDen << " " << maxDen << exit(FatalError);
        consistentRAU = volScalarField::New("consistentRAU", 1.0/denominator);
        const volVectorField oldForce
        (
            fvc::grad(p)
          + fvc::reconstruct(ghf*fvc::snGrad(rho)*mesh.magSf())
        );
        const bool audit = simple.dict().getOrDefault<bool>("consistentAudit", false);
        tmp<volVectorField> originalVelocity(nullptr);
        if (audit) originalVelocity = HbyA - rAU*oldForce;
        HbyA = constrainHbyA
        (
            HbyA + (consistentRAU()-rAU)*oldForce, U, p
        );
        // Use the actual pressure matrix flux, including its non-orthogonal
        // and boundary terms. A differently approximated face gradient would
        // not establish the same discrete fixed point. Coefficient interpolation
        // must be linear; consistentAudit checks both identities when enabled.
        fvScalarMatrix deltaPressureEqn
        (
            fvm::laplacian(rho*(consistentRAU()-rAU), p)
        );
        tmp<surfaceScalarField> deltaPressureFlux(deltaPressureEqn.flux());
        if (audit)
        {
            const volVectorField samePressureVelocity(HbyA-consistentRAU()*oldForce);
            const scalar velocityDefect = gMax
            (
                mag(samePressureVelocity.primitiveField()
                  - originalVelocity().primitiveField())
            )/Foam::max(gMax(mag(originalVelocity().primitiveField())), VSMALL);
            fvScalarMatrix oldPressureEqn(fvm::laplacian(rho*rAU, p));
            fvScalarMatrix newPressureEqn(fvm::laplacian(rho*consistentRAU(), p));
            const surfaceScalarField oldPressureFlux(oldPressureEqn.flux());
            const surfaceScalarField newPressureFlux(newPressureEqn.flux());
            const surfaceScalarField fluxDefect
            (
                newPressureFlux-oldPressureFlux-deltaPressureFlux()
            );
            scalar fluxError = gMax(mag(fluxDefect.primitiveField()));
            scalar fluxScale = Foam::max(gMax(mag(newPressureFlux.primitiveField())), VSMALL);
            forAll(fluxDefect.boundaryField(), patchi)
            {
                if (fluxDefect.boundaryField()[patchi].size())
                {
                    fluxError = Foam::max(fluxError, Foam::max(mag(fluxDefect.boundaryField()[patchi])));
                    fluxScale = Foam::max(fluxScale, Foam::max(mag(newPressureFlux.boundaryField()[patchi])));
                }
            }
            reduce(fluxError, maxOp<scalar>());
            reduce(fluxScale, maxOp<scalar>());
            const scalar fluxRelativeDefect = fluxError/fluxScale;
            Info<< "MEMBRANE_CONSISTENT_AUDIT velocity_relative=" << velocityDefect
                << " flux_relative=" << fluxRelativeDefect << endl;
            if (!(velocityDefect <= 1e-10) || !(fluxRelativeDefect <= 1e-10))
                FatalErrorInFunction << "Consistent discrete fixed-point audit failed"
                    << exit(FatalError);
        }
        phi += deltaPressureFlux();
    }
    const volScalarField& pressureRAU = simple.consistent() ? consistentRAU() : rAU;

    // Solve the equation
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*pressureRAU, p) == fvc::div(phi)
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
    U = HbyA - pressureRAU*
    (
        fvc::grad(p)
      + fvc::reconstruct(ghf*fvc::snGrad(rho)*mesh.magSf())
    );
    U.correctBoundaryConditions();

}
