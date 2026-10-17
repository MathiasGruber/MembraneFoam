rho = rho0 * (1.0 + rho_mACoeff * m_A);

volScalarField rhoD_AB
(
    "rho*D_AB",
    rho * max(D_AB_Coeff * (1.0 - D_AB_mACoeff * m_A), D_AB_Min)
);

fvScalarMatrix m_AEqn
(
    fvm::div(phi, m_A)
  - fvm::laplacian(rhoD_AB, m_A)
);

// Relax to ensure diagonal-dominance
m_AEqn.relax();

// Solve for solute
solve(m_AEqn);
