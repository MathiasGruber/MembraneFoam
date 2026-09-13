"""2012 chamber assembly and shared chamber operating conditions."""

import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from .channel import generate, header

ROOT = Path(__file__).resolve().parents[2]


def initialize_flow(case, *, smoother=None):
    """Compute a flow guess without over-solving the initial pressure equation."""
    case = Path(case)
    solution = case / "system/fvSolution"
    original = solution.read_text()
    try:
        settings = original.replace("tolerance 1e-12", "tolerance 1e-10")
        if smoother is not None:
            settings = re.sub(r"\bsmoother\s+\w+;", f"smoother {smoother};", settings)
        solution.write_text(settings)
        with (case / "log.potentialSalt").open("w") as log:
            subprocess.run(
                ["potentialSalt", "-case", str(case)],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    finally:
        solution.write_text(original)


def configure(case, chamber, flow=50, draw=1, end=24000):
    case = Path(case).resolve()
    if "Mesh OK" not in (case / "log.checkMesh").read_text():
        raise ValueError("A passing checkMesh is required")
    if (case / "0").exists():
        raise FileExistsError("0 exists; use a fresh case")
    (case / "0").mkdir()

    def put(name, text):
        (case / name).write_text(text)

    sets = []
    patches = []
    for name, ylo, yhi, zlo, zhi in [
        ("feedIn", -1, 0, -1, 0),
        ("feedOut", 0, 1, -1, 0),
        ("drawIn", 0, 1, 0, 1),
        ("drawOut", -1, 0, 0, 1),
    ]:
        sets += [
            f"{{ name {name}; type faceSet; action new; source patchToFace; patches (inlet); }}",
            f"{{ name {name}; type faceSet; action subset; source boxToFace; box (-1 {ylo} {zlo}) (1 {yhi} {zhi}); }}",
        ]
        patches += [
            f"{{ name {name}; patchInfo {{ type patch; }} constructFrom set; set {name}; }}"
        ]
    put(
        "system/topoSetDict",
        header("topoSetDict") + "actions (\n" + "\n".join(sets) + "\n);\n",
    )
    put(
        "system/createPatchDict",
        header("createPatchDict")
        + "pointSync false; patches (\n"
        + "\n".join(patches)
        + "\n);\n",
    )
    for tool in ["topoSet", "createPatch"]:
        with (case / ("log." + tool + "-ports")).open("w") as log:
            subprocess.run(
                [tool, "-case", str(case)]
                + (["-overwrite"] if tool == "createPatch" else []),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    # Half-width domain: half of the published volumetric flow per compartment.
    q = flow * 1e-6 / 120
    c = draw * 58.44
    mf = 2 * c / (997.1 + math.sqrt(997.1**2 + 4 * 694 * c))
    put(
        "0/m_A",
        header("m_A", "volScalarField")
        + f"""
    dimensions [0 0 0 0 0 0 0]; internalField uniform 0;
    boundaryField {{
     feedIn {{ type fixedValue; value uniform 0; }}
     drawIn {{ type fixedValue; value uniform {mf}; }}
     membrane {{ type explicitFOmembraneSolute; value uniform 0; }}
     Symmetri {{ type symmetryPlane; }}
     ".*" {{ type zeroGradient; }}
    }}
    """,
    )
    put(
        "0/U",
        header("U", "volVectorField")
        + f"""
    dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0);
    boundaryField {{
     "(feedIn|drawIn)" {{ type flowRateInletVelocity; volumetricFlowRate {q}; value uniform (0 0 0); }}
     "(feedOut|drawOut)" {{ type zeroGradient; }}
     membrane {{ type explicitFOmembraneVelocity; forwardDirection (0 0 1); eq advanced; value uniform (0 0 0); }}
     Symmetri {{ type symmetryPlane; }}
     "(fixedWalls|cylinderInlets|carvingTop|inlet_connector)" {{ type fixedValue; value uniform (0 0 0); }}
    }}
    """,
    )
    put(
        "0/p",
        header("p", "volScalarField")
        + """
    dimensions [1 -1 -2 0 0 0 0]; internalField uniform 0;
    boundaryField {
     "(feedOut|drawOut)" { type fixedValue; value uniform 0; }
     Symmetri { type symmetryPlane; }
     "(feedIn|drawIn|membrane|fixedWalls|cylinderInlets|carvingTop|inlet_connector)" { type zeroGradient; }
    }
    """,
    )
    put(
        "system/setFieldsDict",
        header("setFieldsDict")
        + f"""
    defaultFieldValues (volScalarFieldValue m_A 0);
    regions (boxToCell {{ box (-1 -1 0) (1 1 1); fieldValues (volScalarFieldValue m_A {mf}); }});
    """,
    )
    # Transport settings are the 2012 values from the generated channel template.
    with tempfile.TemporaryDirectory() as tmp:
        template = Path(tmp) / "template"
        generate(template, flow=flow, draw=draw, end=end)
        for name in [
            "constant/transportProperties",
            "constant/g",
            "system/controlDict",
            "system/fvSchemes",
            "system/fvSolution",
        ]:
            shutil.copyfile(template / name, case / name)
    # Non-orthogonal corrections for the curved-port chamber mesh.
    s = (
        (case / "system/fvSolution")
        .read_text()
        .replace("nNonOrthogonalCorrectors 0", "nNonOrthogonalCorrectors 2")
        .replace("tolerance 1e-11;", "tolerance 1e-12;")
        .replace("m_A 1e-7;", "m_A 1e-9;")
        .replace("m_A 0.7;", "m_A 0.2;")
        .replace("smoother DICGaussSeidel;", "smoother GaussSeidel;")
    )
    if chamber == "a":
        s = s.replace("p 0.3;", "p 0.15;").replace("U 0.7;", "U 0.3;")
        s = s.replace(
            "smoother GaussSeidel;",
            "smoother GaussSeidel; nCellsInCoarsestLevel 100; directSolveCoarsest yes;",
        )
    put("system/fvSolution", s)
    schemes = (
        (case / "system/fvSchemes")
        .read_text()
        .replace("Gauss limitedLinear 1", "Gauss vanLeer")
    )
    put("system/fvSchemes", schemes)
    put(
        "case.json",
        json.dumps(
            dict(
                paper="10.3390/membranes2040764",
                geometry="2012 chamber",
                chamber=chamber.upper(),
                flow_ml_min=flow,
                draw_molar=draw,
                half_width=True,
            ),
            indent=2,
        )
        + "\n",
    )
    for tool in ["setFields"]:
        with (case / ("log." + tool)).open("w") as log:
            subprocess.run(
                [tool, "-case", str(case)],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    initialize_flow(case, smoother="DICGaussSeidel" if chamber == "a" else None)


def prepare(case, chamber, end=24000, refinement=1):
    """Generate, assemble, check, and initialize a chamber in a fresh directory."""
    if not math.isfinite(refinement) or refinement < 1:
        raise ValueError("Refinement must be finite and >= 1")
    case = Path(case).resolve()
    (case / "system").mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory() as tmp:
        template = Path(tmp) / "case"
        generate(template, end=end)
        for name in ("controlDict", "fvSchemes", "fvSolution"):
            shutil.copyfile(template / "system" / name, case / "system" / name)
    with (case / "log.generator").open("w") as log:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / f"examples/papers/2012/generators/chamber_{chamber}.py"),
            ],
            cwd=case,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )

    if refinement != 1:
        mesh = case / "system/blockMeshDict"
        pattern = r"(\bhex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)"

        def refined(match):
            counts = [math.ceil(int(match[i]) * refinement) for i in (2, 3, 4)]
            return match[1] + "(" + " ".join(map(str, counts)) + ")"

        text, count = re.subn(pattern, refined, mesh.read_text())
        if not count:
            raise ValueError("No block cell counts found")
        mesh.write_text(text)

    def tool(name, *args, append=False):
        with (case / ("log." + name)).open("a" if append else "w") as log:
            subprocess.run(
                [name, *args],
                cwd=case,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    tool("blockMesh")
    if refinement != 1:
        tool("checkMesh")
        if (case / "constant/polyMesh/sets/skewFaces").exists():
            # The legacy port intersections can create sliver faces as counts
            # change. Repair the quarter mesh before mirroring to keep membrane
            # partners identical. Retain the final, unmodified checkMesh gate.
            template = (
                Path(os.environ["WM_PROJECT_DIR"])
                / "etc/caseDicts/annotated/collapseDict"
            )
            controls = template.read_text()
            controls = re.sub(
                r"minimumEdgeLength\s+[^;]+;", "minimumEdgeLength 1e-8;", controls
            )
            controls = re.sub(
                r"maximumMergeAngle\s+[^;]+;", "maximumMergeAngle 1;", controls
            )
            controls = re.sub(
                r'#include\s+"meshQualityDict";',
                '#includeEtc "caseDicts/meshQualityDict"\n minVol 1e-18;',
                controls,
            )
            (case / "system/collapseDict").write_text(controls)
            tool("collapseEdges", "-collapseFaceSet", "skewFaces", "-overwrite")
            # Remove mesh-filter diagnostic fields before physical initialization.
            shutil.rmtree(case / "0")
    for normal in ("0 0 1", "0 1 0"):
        (case / "system/mirrorMeshDict").write_text(
            header("mirrorMeshDict")
            + f"""
planeType pointAndNormal;
pointAndNormalDict {{ basePoint (0 0 0); normalVector ({normal}); }}
planeTolerance 1e-8;
"""
        )
        tool("mirrorMesh", "-overwrite", append=True)
    (case / "system/topoSetDict").write_text(
        header("topoSetDict")
        + """
actions (
 { name membraneFaces; type faceSet; action new; source boxToFace; box (-1 -1 -1e-10) (1 1 1e-10); }
 { name membraneZone; type faceZoneSet; action new; source setToFaceZone; faceSet membraneFaces; }
);
"""
    )
    tool("topoSet")
    (case / "system/createBafflesDict").write_text(
        header("createBafflesDict")
        + """
internalFacesOnly true;
baffles {
 membrane { type faceZone; zoneName membraneZone;
 patches { master { name membrane; type wall; } slave { name membrane; type wall; } }
 }
}
"""
    )
    tool("createBaffles", "-overwrite")
    tool("checkMesh")
    configure(case, chamber, end=end)
    metadata = case / "case.json"
    record = json.loads(metadata.read_text())
    record["refinement"] = refinement
    metadata.write_text(json.dumps(record, indent=2) + "\n")



def configure_2016(path, flow, draw, feed, resistance, end):
    """Apply the 2016 FO parameters to a generated case with standard patch names."""
    properties = (path / "constant/transportProperties").read_text()
    for key, value in [("A", 1.61111e-12), ("B", 8.33333e-8), ("K", resistance)]:
        properties = re.sub(rf"^{key} .*;", f"{key} {value};", properties, flags=re.M)
    (path / "constant/transportProperties").write_text(properties)
    control = path / "system/controlDict"
    control.write_text(
        control.read_text()
        .replace(f"writeInterval {end}", "writeInterval 1000")
        .replace("writeFormat ascii", "writeFormat binary")
    )
    # Prescribed volumetric flow is half the whole-module experimental flow.
    q = flow * 1e-6 / 60 / 2
    for field in ("U", "p", "m_A"):
        p = path / "0" / field
        s = p.read_text().replace("type empty;", "type symmetryPlane;")
        initial = "(0 0 0)" if field == "U" else str(feed if field == "m_A" else 0)
        s = re.sub(
            r"internalField\s+nonuniform\s+List<[^>]+>\s+\d+\s*\(.*?\);",
            f"internalField uniform {initial};",
            s,
            flags=re.S,
        )
        if field == "U":
            for name in ("feedIn", "drawIn"):
                s = re.sub(
                    rf"{name}\s*\{{[^}}]*\}}",
                    f"{name} {{ type flowRateInletVelocity; volumetricFlowRate constant {q}; value uniform (0 0 0); }}",
                    s,
                )
        elif field == "m_A":
            for name, value in [("feedIn", feed), ("drawIn", draw)]:
                s = re.sub(
                    rf"{name}\s*\{{[^}}]*\}}",
                    f"{name} {{ type fixedValue; value uniform {value}; }}",
                    s,
                )
        p.write_text(s)
    (path / "system/setFieldsDict").write_text(
        header("setFieldsDict")
        + f"""
defaultFieldValues (volScalarFieldValue m_A {feed});
regions (boxToCell {{ box (-1 -1 0) (1 1 1); fieldValues (volScalarFieldValue m_A {draw}); }});
"""
    )
    p = path / "system/fvSolution"
    s = (
        p.read_text()
        .replace("nNonOrthogonalCorrectors 0", "nNonOrthogonalCorrectors 2")
        .replace("m_A 1e-7", "m_A 1e-9")
        .replace("SIMPLE {", "SIMPLE { maxRelativeSaltImbalance 0.001;")
        .replace("tolerance 1e-11", "tolerance 1e-12")
        .replace("m_A 0.7", "m_A 0.3")
    )
    p.write_text(s)
    p = path / "system/fvSchemes"
    p.write_text(p.read_text().replace("Gauss limitedLinear 1", "Gauss vanLeer"))
