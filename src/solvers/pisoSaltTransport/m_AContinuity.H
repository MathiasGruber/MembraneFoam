{
    scalar m_A_influx_old = m_A_influx;
    scalar m_A_mass_old = m_A_mass;

    m_A_mass = fvc::domainIntegrate(m_A * rho0 * rho_mACoeff).value();
    Info << m_A.name() << " mass: old = " << m_A_mass_old << "    new = " << m_A_mass << endl;
    m_A_influx = 0.0;

    forAll(m_A.boundaryField(), patchi)
    {
        m_A_influx -= sum(phi.boundaryField()[patchi] / rho.boundaryField()[patchi] * m_A.boundaryField()[patchi]);
    }
    m_A_influx *= rho0.value() * rho_mACoeff.value();
    Info << m_A.name() << " influx: old = " << m_A_influx_old << "    new = " << m_A_influx << endl;
}
