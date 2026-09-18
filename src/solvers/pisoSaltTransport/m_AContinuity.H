// SPDX-License-Identifier: GPL-3.0-or-later
{
    // Match Time::run(): accumulated time may lie just below the endpoint.
    if
    (
        runTime.writeTime()
     || runTime.value() >= runTime.endTime().value() - 0.5*runTime.deltaTValue()
    )
    {
        const scalar saltMass = fvc::domainIntegrate(rho*m_A).value();
        const scalar saltStorage = fvc::domainIntegrate(fvc::ddt(rho, m_A)).value();
        const scalar massStorage = fvc::domainIntegrate(fvc::ddt(rho)).value();
        const auto saltFlux = membraneSaltFlux
        (
            mesh, phi, m_A, rho, D_AB_Coeff, D_AB_mACoeff, D_AB_Min
        );
        scalar massOut = 0, saltOut = 0, throughput = 0, membraneTransfer = 0;
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled()) continue;
            massOut += gSum(phi.boundaryField()[patchi]);
            throughput += gSum(mag(phi.boundaryField()[patchi]))/2;
            saltOut += gSum(saltFlux().boundaryField()[patchi]);
            const word type = U.boundaryField()[patchi].type();
            if (type == "explicitFOmembraneVelocity" || type == "explicitROmembraneVelocity")
                membraneTransfer += gSum(mag(saltFlux().boundaryField()[patchi]))/2;
        }
        Info<< "MEMBRANE_TRANSIENT time=" << runTime.value()
            << " salt_mass_kg=" << saltMass
            << " salt_storage_kg_s=" << saltStorage
            << " salt_outflow_kg_s=" << saltOut
            << " salt_balance_kg_s=" << saltStorage+saltOut
            << " relative_salt_balance=" << mag(saltStorage+saltOut)/Foam::max(membraneTransfer, VSMALL)
            << " mass_storage_kg_s=" << massStorage
            << " mass_outflow_kg_s=" << massOut
            << " relative_mass_balance=" << mag(massStorage+massOut)/Foam::max(throughput, VSMALL)
            << nl;
    }
}
