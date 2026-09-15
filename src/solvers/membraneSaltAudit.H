// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef membraneSaltAudit_H
#define membraneSaltAudit_H
#include "fvCFD.H"
#include <algorithm>

namespace Foam
{
inline tmp<surfaceScalarField> membraneSaltFlux
(
    fvMesh& mesh, const surfaceScalarField& phi, const volScalarField& m,
    const volScalarField& rho, const dimensionedScalar& diffusivity,
    const dimensionedScalar& concentrationCoefficient,
    const dimensionedScalar& minimumDiffusivity
)
{
    mesh.setFluxRequired(m.name());
    const volScalarField rhoD
    (
        "auditRhoD", rho*max(diffusivity*(1-concentrationCoefficient*m), minimumDiffusivity)
    );
    fvScalarMatrix audit(fvm::div(phi, m) - fvm::laplacian(rhoD, m));
    return audit.flux();
}

inline scalar relativeMembraneSaltImbalance
(
    const surfaceScalarField& saltFlux, const volVectorField& U
)
{
    scalar total = 0, transfer = 0;
    forAll(U.mesh().boundary(), patchi)
    {
        if (U.mesh().boundary()[patchi].coupled()) continue;
        total += gSum(saltFlux.boundaryField()[patchi]);
        const word type = U.boundaryField()[patchi].type();
        if (type == "explicitFOmembraneVelocity" || type == "explicitROmembraneVelocity")
            transfer += gSum(mag(saltFlux.boundaryField()[patchi]))/2;
    }
    return mag(total)/std::max<scalar>(transfer, VSMALL);
}
}
#endif
