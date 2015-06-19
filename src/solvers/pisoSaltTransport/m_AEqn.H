
volScalarField rhoD_AB
(
    "rho*D_AB",
    rho * max(D_AB_Coeff * (1.0 - D_AB_mACoeff * m_A), D_AB_Min)
);

solve
(
    fvm::ddt(rho,m_A)
  + fvm::div(phi, m_A)
  - fvm::laplacian(rhoD_AB, m_A)
); 

Info << m_A.name() <<  ": min = " << min(m_A).value() << "    max = " << max(m_A).value() << endl;
