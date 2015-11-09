rho = rho0 * (1.0 + rho_mACoeff * m_A);

// Calculate rho * D_AB
volScalarField rhoD_AB
(
    "rho*D_AB*Diff_ratio",
    rho * max(D_AB_Coeff * (1.0 - D_AB_mACoeff * m_A), D_AB_Min) * Diff_ratio
);

// Calculate phiU inside vesicle
surfaceScalarField phiUvesicle
(
    "phiUvesicle",
    phi * fvc::interpolate(inVesicle)
);

// Solute equation
fvScalarMatrix m_AEqn
(
    fvm::div(phiUvesicle, m_A)
  - fvm::laplacian(rhoD_AB, m_A)
);

// Relax to ensure diagonal-dominance
m_AEqn.relax();

// Solve for solute
solve(m_AEqn);