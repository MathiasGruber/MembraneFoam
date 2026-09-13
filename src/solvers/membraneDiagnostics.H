// Integral checks for the converged iterate (SI units). Physical boundary
// transport includes convection and diffusion; processor boundaries are skipped.
{
    // Audit the actual finite-volume salt flux, including non-orthogonal
    // corrections, separately from the physical boundary-gradient estimate.
    const surfaceScalarField numericalSaltFlux
    (
        membraneSaltFlux(mesh, phi, m_A, rho, D_AB_Coeff, D_AB_mACoeff, D_AB_Min)
    );
    scalar numericalSaltBalance = 0, membraneSaltTransfer = 0;
    scalar totalMass = 0, totalSalt = 0, throughput = 0;
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].coupled()) continue;
        const auto& massFlux = phi.boundaryField()[patchi];
        const auto& solute = m_A.boundaryField()[patchi];
        const auto& area = mesh.boundary()[patchi].magSf();
        const scalarField diffusivity
        (
            rho.boundaryField()[patchi]
           *max(D_AB_Coeff.value()*(1-D_AB_mACoeff.value()*solute),D_AB_Min.value())
        );
        const scalar mass = gSum(massFlux);
        const scalar salt = gSum(massFlux*solute - diffusivity*solute.snGrad()*area);
        const scalar numericalSalt = gSum(numericalSaltFlux.boundaryField()[patchi]);
        numericalSaltBalance += numericalSalt;
        Info<< "MEMBRANE_DISCRETE patch=" << mesh.boundary()[patchi].name()
            << " salt_kg_s=" << numericalSalt << nl;
        const word velocityType = U.boundaryField()[patchi].type();
        if (velocityType == "explicitFOmembraneVelocity"
            || velocityType == "explicitROmembraneVelocity")
        {
            membraneSaltTransfer += gSum(mag(numericalSaltFlux.boundaryField()[patchi]))/2;
            const scalar doubleArea = std::max<scalar>(gSum(area), VSMALL);
            Info<< "MEMBRANE_FLUX patch=" << mesh.boundary()[patchi].name()
                << " water_kg_m2_h=" << gSum(mag(massFlux))*3600/doubleArea
                << " salt_g_m2_h=" << gSum(mag(massFlux*solute-diffusivity*solute.snGrad()*area))*3.6e6/doubleArea << nl;
            OFstream surface(runTime.path()/(mesh.boundary()[patchi].name()+".csv"));
            surface.precision(12);
            surface << "x_m,y_m,z_m,nx,ny,nz,mass_fraction,normal_velocity_m_s,area_m2" << nl;
            const auto normals = mesh.boundary()[patchi].nf();
            const auto& centres = mesh.boundary()[patchi].Cf();
            forAll(solute,i)
                surface << centres[i].x() << ',' << centres[i].y() << ',' << centres[i].z()
                    << ',' << normals()[i].x() << ',' << normals()[i].y() << ',' << normals()[i].z()
                    << ',' << solute[i] << ',' << (U.boundaryField()[patchi][i] & normals()[i]) << "," << area[i] << nl;
        }
        totalMass += mass;
        totalSalt += salt;
        throughput += gSum(mag(massFlux));
        Info<< "MEMBRANE_CHECK patch=" << mesh.boundary()[patchi].name()
            << " mass_kg_s=" << mass << " salt_kg_s=" << salt << nl;
    }
    Info<< "MEMBRANE_CHECK total_mass_kg_s=" << totalMass
        << " total_salt_kg_s=" << totalSalt
        << " relative_mass_imbalance=" << mag(totalMass)/std::max<scalar>(throughput/2, VSMALL) << nl;
    Info<< "MEMBRANE_DISCRETE total_salt_kg_s=" << numericalSaltBalance
        << " membrane_salt_transfer_kg_s=" << membraneSaltTransfer
        << " relative_salt_imbalance=" << mag(numericalSaltBalance)/std::max<scalar>(membraneSaltTransfer, VSMALL) << nl;
}
