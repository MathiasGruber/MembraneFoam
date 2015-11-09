// Output some information for the user
Info<< "rho max/min : "
    << max(rho).value() << " / "
    << min(rho).value() << endl;

Info<< "pressure max/min : "
    << max(p).value() << " / "
    << min(p).value() << endl;

Info<< "m_A max/min : "
    << max(m_A).value() << " / "
    << min(m_A).value() << endl;

Info<< "U max/min : "
    << max(U).value() << " / "
    << min(U).value() << endl;

Info<< "mu*laplacian(U) max/min : "
    << max(mu_laplacianU).value() << " / "
    << min(mu_laplacianU).value() << endl;

// Calculate new density
rho = rho0 * (1.0 + rho_mACoeff * m_A);

volVectorField HbyA("HbyA", U);

if (pressureImplicitPorosity)
{
    HbyA = trTU() & UEqn().H();
}
else
{
    HbyA = trAU()*UEqn().H();
}

UEqn.clear();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho*HbyA) & mesh.Sf()
);

fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn =
        (
            fvm::laplacian(rho*trTU(), p)
         ==
            fvc::div(phiHbyA)
        );
    }
    else
    {
        tpEqn =
        (
            fvm::laplacian(rho*trAU(), p)
         ==
            fvc::div(phiHbyA)
        );
    }

    tpEqn().setReference(pRefCell, pRefValue);

    fvOptions.constrain(tpEqn(), rho.name());

    tpEqn().solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - tpEqn().flux();
    }
}

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U = HbyA - (trTU() & fvc::grad(p));
}
else
{
    U = HbyA - trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();
fvOptions.correct(U);