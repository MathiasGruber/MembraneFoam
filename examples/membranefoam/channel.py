#!/usr/bin/env python3
"""Generate a conformal two-sided membrane channel; no legacy mesh tools needed.

The 2012 preset uses published chamber-A compartment dimensions and membrane
properties. Full-width uniform slots replace its shaped ports, so this is a
controlled channel benchmark, not a claim of reproducing chamber-A figures.
"""

from pathlib import Path


def header(name, cls="dictionary"):
    return f"FoamFile {{ version 2.0; format ascii; class {cls}; object {name}; }}\n"


def generate(path, nx=60, nz=24, flow=20.0, draw=0.5, end=12000, mode="FO"):
    if nx < 1 or nz < 2 or end < 1:
        raise ValueError("nx >= 1, nz >= 2 and end >= 1 are required")
    if not (flow > 0 and 0 <= draw <= 1.5):
        raise ValueError(
            "Use positive flow and draw between 0 and 1.5 M (property-fit range)"
        )
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"{path} exists; choose a new run directory")
    for d in ("0", "system", "constant"):
        (path / d).mkdir(parents=True, exist_ok=True)

    def put(name, text):
        (path / name).write_text(text)

    # SI: 0.44 LMH/bar, 0.087 LMH, 0.72 s/um (2012 section 3.1).
    a, b, k = 0.44 / 3.6e11, 0.087 / 3.6e6, 0.72e6
    length, width, height = 0.03, 0.015, 0.001
    speed = flow * 1e-6 / 60 / (width * height)
    # Convert mol/L NaCl to mass fraction using rho=997.1+694*m.
    concentration = draw * 58.44
    mf = 2 * concentration / (997.1 + (997.1**2 + 4 * 694 * concentration) ** 0.5)
    points = []
    for z0, z1 in ((-height, 0), (0, height)):
        points += [
            (0, 0, z0),
            (length, 0, z0),
            (length, width, z0),
            (0, width, z0),
            (0, 0, z1),
            (length, 0, z1),
            (length, width, z1),
            (0, width, z1),
        ]
    point_text = "\n".join(f"({x} {y} {z})" for x, y, z in points)
    put(
        "system/blockMeshDict",
        header("blockMeshDict")
        + f"""
scale 1;
mergeType topology;
vertices ({point_text});
blocks (
 hex (0 1 2 3 4 5 6 7) ({nx} 1 {nz}) simpleGrading (1 1 0.1)
 hex (8 9 10 11 12 13 14 15) ({nx} 1 {nz}) simpleGrading (1 1 10)
);
edges ();
boundary (
 feedIn {{ type patch; faces ((0 4 7 3)); }}
 feedOut {{ type patch; faces ((1 2 6 5)); }}
 drawIn {{ type patch; faces ((9 10 14 13)); }}
 drawOut {{ type patch; faces ((8 12 15 11)); }}
 membrane {{ type wall; faces ((4 5 6 7) (8 11 10 9)); }}
 walls {{ type wall; faces ((0 3 2 1) (12 13 14 15)); }}
 sides {{ type empty; faces ((0 1 5 4) (3 7 6 2) (8 9 13 12) (11 15 14 10)); }}
);
mergePatchPairs ();
""",
    )
    put(
        "constant/transportProperties",
        header("transportProperties")
        + f"""
pi_mACoeff pi_mACoeff [1 -1 -2 0 0 0 0] 80510000;
mu0 mu0 [1 -1 -1 0 0 0 0] 0.00089;
mu_mACoeff mu_mACoeff [0 0 0 0 0 0 0] 1.63;
D_AB_Min D_AB_Min [0 2 -1 0 0 0 0] 1.45e-9;
D_AB_Coeff D_AB_Coeff [0 2 -1 0 0 0 0] 1.61e-9;
D_AB_mACoeff D_AB_mACoeff [0 0 0 0 0 0 0] 14;
rho0 rho0 [1 -3 0 0 0 0 0] 997.1;
rho_mACoeff rho_mACoeff [0 0 0 0 0 0 0] {694 / 997.1};
A {a};
B {b};
K {k if mode == "FO" else a};
""",
    )
    put(
        "constant/g",
        header("g", "uniformDimensionedVectorField")
        + "dimensions [0 1 -2 0 0 0 0]; value (0 0 0);\n",
    )
    solute = "explicitFOmembraneSolute" if mode == "FO" else "explicitROmembraneSolute"
    velocity = (
        "explicitFOmembraneVelocity" if mode == "FO" else "explicitROmembraneVelocity"
    )
    feed_m, draw_m = (0.0, mf) if mode == "FO" else (mf, 0.0)
    values = "\n".join([str(feed_m)] * (nx * nz) + [str(draw_m)] * (nx * nz))
    put(
        "0/m_A",
        header("m_A", "volScalarField")
        + f"""
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar> {2 * nx * nz} ({values});
boundaryField {{
 feedIn {{ type fixedValue; value uniform {feed_m}; }}
 drawIn {{ type fixedValue; value uniform {draw_m}; }}
 "(feedOut|drawOut)" {{ type zeroGradient; }}
 membrane {{ type {solute}; R 0.99; value uniform {mf / 2}; }}
 walls {{ type zeroGradient; }}
 sides {{ type empty; }}
}}
""",
    )
    uvalues = "\n".join(
        [f"({speed} 0 0)"] * (nx * nz) + [f"({-speed} 0 0)"] * (nx * nz)
    )
    put(
        "0/U",
        header("U", "volVectorField")
        + f"""
dimensions [0 1 -1 0 0 0 0]; internalField nonuniform List<vector> {2 * nx * nz} ({uvalues});
boundaryField {{
 feedIn {{ type fixedValue; value uniform ({speed} 0 0); }}
 drawIn {{ type fixedValue; value uniform ({-speed} 0 0); }}
 "(feedOut|drawOut)" {{ type zeroGradient; }}
 membrane {{ type {velocity}; K {a}; forwardDirection (0 0 1); eq advanced; value uniform (0 0 0); }}
 walls {{ type fixedValue; value uniform (0 0 0); }}
 sides {{ type empty; }}
}}
""",
    )
    pvalues = "\n".join(
        [str(6e6 if mode == "RO" else 0)] * (nx * nz) + ["0"] * (nx * nz)
    )
    put(
        "0/p",
        header("p", "volScalarField")
        + f"""
dimensions [1 -1 -2 0 0 0 0]; internalField nonuniform List<scalar> {2 * nx * nz} ({pvalues});
boundaryField {{
 "(feedIn|drawIn|membrane|walls)" {{ type zeroGradient; }}
 feedOut {{ type fixedValue; value uniform {6e6 if mode == "RO" else 0}; }}
 drawOut {{ type fixedValue; value uniform 0; }}
 sides {{ type empty; }}
}}
""",
    )
    put(
        "system/controlDict",
        header("controlDict")
        + f"""
application simpleSaltTransport;
startFrom startTime; startTime 0; stopAt endTime; endTime {end}; deltaT 1;
writeControl timeStep; writeInterval {end}; purgeWrite 2;
writeFormat ascii; writePrecision 12; runTimeModifiable false;
libs ("libDHIBoundaryConditions.so");
""",
    )
    put(
        "system/fvSchemes",
        header("fvSchemes")
        + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes {
 default none;
 div(phi,U) Gauss linearUpwind grad(U);
 div(phi,m_A) Gauss limitedLinear 1;
 div((mu*dev2(T(grad(U))))) Gauss linear;
 div((mu*dev2(grad(U).T()))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
""",
    )
    put(
        "system/fvSolution",
        header("fvSolution")
        + """
solvers {
 p { solver GAMG; tolerance 1e-11; relTol 0.001; smoother DICGaussSeidel; }
 pFinal { $p; relTol 0; }
 "(U|UFinal|m_A)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-15; relTol 0; }
}
SIMPLE {
 nNonOrthogonalCorrectors 0;
 residualControl { p 1e-6; U 1e-6; m_A 1e-7; }
}
PISO { nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0; momentumPredictor true; }
potentialFlow { nNonOrthogonalCorrectors 10; }
relaxationFactors { fields { p 0.3; } equations { U 0.7; m_A 0.7; } }
""",
    )
    put(
        "system/decomposeParDict",
        header("decomposeParDict")
        + """
numberOfSubdomains 4;
method simple;
simpleCoeffs { n (4 1 1); delta 0.001; }
""",
    )
    put(
        "case.json",
        __import__("json").dumps(
            dict(
                mode=mode,
                nx=nx,
                nz=nz,
                flow_ml_min=flow,
                draw_molar=draw,
                paper="10.3390/membranes2040764",
                geometry="idealized uniform-slot channel",
                membrane_area_m2=length * width,
            ),
            indent=2,
        )
        + "\n",
    )
