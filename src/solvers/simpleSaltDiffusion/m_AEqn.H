rho = rho0 * (1.0 + rho_mACoeff * m_A);

// Calculate rho * D_AB
volScalarField rhoD_AB
(
    "rho*D_AB*Diff_ratio",
    rho * max(D_AB_Coeff * (1.0 - D_AB_mACoeff * m_A), D_AB_Min) * Diff_ratio
);

// Solute equation
fvScalarMatrix m_AEqn
(
  - fvm::laplacian(rhoD_AB, m_A)
);

// Relax to ensure diagonal-dominance
m_AEqn.relax();

// Solve for solute
solve(m_AEqn);

// Output info
Info<< "m_A max/min : "
    << max(m_A).value() << " / "
    << min(m_A).value() << endl;