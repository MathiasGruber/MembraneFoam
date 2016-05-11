// Pre-OF3: #include "readPISOControls.H"
// Now directly in this file:
const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
const int nOuterCorr =
    pisoDict.lookupOrDefault<int>("nOuterCorrectors", 1);
const int nCorr =
    pisoDict.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr =
    pisoDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
const bool momentumPredictor =
    pisoDict.lookupOrDefault("momentumPredictor", true);
const bool transonic =
    pisoDict.lookupOrDefault("transonic", false);

// Get the time controls
#include "readTimeControls.H"
