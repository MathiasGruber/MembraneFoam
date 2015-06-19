rho = rho0 * (1.0 + rho_mACoeff * m_A);

{
    // OpenFoam 2.2 and below:
    // volScalarField rUA = 1.0/UEqn.A();
    // surfaceScalarField rUAf = fvc::interpolate(rUA);
    tmp<volScalarField> temp0 = 1.0/UEqn.A();
    volScalarField rUA = temp0();
    
    tmp<surfaceScalarField> temp1 = fvc::interpolate(rUA);
    surfaceScalarField rUAf = temp1();

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(rho) *
        (
            (fvc::interpolate(U) & mesh.Sf())
          + rUAf*fvc::ddtCorr(rho, U, phi)
        )
    );

    phi = phiU - fvc::interpolate(rho)*ghf*fvc::snGrad(rho)*rUAf*mesh.magSf();

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rUA, p) 
         == fvc::div(phi) 
          + fvc::ddt(rho)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve
        (
            mesh.solver
            (
                p.name()
              + ((corr == nCorr-1 && nonOrth == nNonOrthCorr) ? "Final" : "")
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    // the standard continuity errors calculated for compressible flow cannot
    // be used here as they rely on a basicThermo object
//    #include "compressibleContinuityErrs.H"

    U += rUA/rho*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
