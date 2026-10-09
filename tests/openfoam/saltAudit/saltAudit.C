// SPDX-License-Identifier: GPL-3.0-or-later
#include "fvCFD.H"
#include "timeSelector.H"
#include "../../../src/solvers/membraneSaltAudit.H"

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    const instantList times = timeSelector::select0(runTime, args);
    if (times.empty()) FatalErrorInFunction << "No field time" << exit(FatalError);
    runTime.setTime(times.last(), times.size() - 1);
    #include "createMesh.H"

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ), mesh
    );
    volScalarField p
    (
        IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ), mesh
    );
    volScalarField m_A
    (
        IOobject("m_A", runTime.timeName(), mesh, IOobject::MUST_READ), mesh
    );
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::MUST_READ), mesh
    );
    IOdictionary properties
    (
        IOobject("transportProperties", runTime.constant(), mesh, IOobject::MUST_READ)
    );
    const dimensionedScalar rho0("rho0", properties);
    const dimensionedScalar beta("rho_mACoeff", properties);
    const dimensionedScalar diffusivity("D_AB_Coeff", properties);
    const dimensionedScalar slope("D_AB_mACoeff", properties);
    const dimensionedScalar minimum("D_AB_Min", properties);
    const volScalarField rho("rho", rho0*(1 + beta*m_A));
    const volScalarField rhoD_AB
    (
        "rho*D_AB", rho*max(diffusivity*(1 - slope*m_A), minimum)
    );
    mesh.setFluxRequired(m_A.name());
    fvScalarMatrix equation(fvm::div(phi, m_A) - fvm::laplacian(rhoD_AB, m_A));
    const auto expected = equation.flux();
    const auto observed = membraneSaltFlux
    (
        mesh, phi, m_A, rho, diffusivity, slope, minimum
    );
    scalar error = sum(mag(expected().primitiveField() - observed().primitiveField()));
    scalar scale = sum(mag(expected().primitiveField()));
    forAll(mesh.boundary(), patchi)
    {
        error += sum(mag(expected().boundaryField()[patchi] - observed().boundaryField()[patchi]));
        scale += sum(mag(expected().boundaryField()[patchi]));
    }
    // Processor-patch counts differ by rank; use a fixed number of collectives.
    reduce(error, sumOp<scalar>());
    reduce(scale, sumOp<scalar>());
    const scalar relative = error/std::max<scalar>(scale, VSMALL);
    Info << "Salt audit/operator relative difference = " << relative << nl;
    if (relative > 1e-12)
    {
        FatalErrorInFunction << "Salt audit does not match the transport operator"
            << exit(FatalError);
    }
    return 0;
}
