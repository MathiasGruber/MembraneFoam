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

// Pre-OF3: #include "readTimeControls.H"
// Now directly in this file:

// Maximum flow Courant number
scalar maxCo(readScalar(pimple.dict().lookup("maxCo")));

// Maximum time scale
scalar maxDeltaT(pimple.dict().lookupOrDefault<scalar>("maxDeltaT", GREAT));