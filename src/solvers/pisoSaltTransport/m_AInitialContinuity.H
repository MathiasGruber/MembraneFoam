scalar m_A_influx = 0.0;
scalar m_A_mass = fvc::domainIntegrate(rho0 * rho_mACoeff * m_A).value();

Info << "m_A: Initial mass = " << m_A_mass << endl << endl;
