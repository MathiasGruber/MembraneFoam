#!/usr/bin/env python3
"""Generate and run a membrane example, exporting validated results automatically."""

import argparse
import json
import math
import re
import os
import subprocess
from pathlib import Path

from membranefoam.channel import generate, header
from membranefoam.chambers import prepare, initialize_flow
from membranefoam.results import collect


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "example", choices=["channel", "chamber-a", "chamber-b", "cf042", "3dblock"]
    )
    parser.add_argument("case", type=Path, help="New case directory")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--nx", type=int, default=60)
    parser.add_argument("--nz", type=int, default=24)
    parser.add_argument("--mode", choices=["FO", "RO"], default="FO")
    parser.add_argument("--end", type=int)
    parser.add_argument(
        "--transient",
        type=float,
        metavar="SECONDS",
        help="Integrate with PISO to this physical time",
    )
    parser.add_argument(
        "--time-step", type=float, help="Fixed transient time step in seconds"
    )
    parser.add_argument(
        "--outer-correctors", type=int, help="PISO outer correctors (default: 4)"
    )
    parser.add_argument(
        "--non-orthogonal-correctors",
        type=int,
        help="PISO non-orthogonal pressure corrections (default: 2)",
    )
    parser.add_argument(
        "--inlets", type=int, default=3, help="3dBlock inlet count (odd, 1–19)"
    )
    parser.add_argument(
        "--angle", type=float, default=45, help="3dBlock inlet angle, degrees"
    )
    parser.add_argument(
        "--refinement",
        type=float,
        default=1,
        help="2012 chamber or 3dBlock mesh refinement factor",
    )
    parser.add_argument(
        "--spacers",
        type=int,
        default=0,
        help="3dBlock middle cylinders per channel (even, 0–18)",
    )
    parser.add_argument(
        "--spacer-radius",
        type=float,
        help="Cylinder radius in metres; required with --spacers",
    )
    parser.add_argument(
        "--flow",
        type=float,
        default=50,
        help="2016 chamber whole-module flow, mL/min per compartment",
    )
    parser.add_argument("--length", type=float, default=None, help="CF042 length, m")
    parser.add_argument("--width", type=float, default=None, help="CF042 width, m")
    parser.add_argument(
        "--height", type=float, default=None, help="CF042 channel height, m"
    )
    parser.add_argument(
        "--spacing", type=float, default=None, help="CF042 tangential mesh spacing, m"
    )
    parser.add_argument(
        "--layers", type=int, default=None, help="CF042 cells through each thin channel"
    )
    parser.add_argument(
        "--offset",
        type=float,
        default=None,
        help="CF042 distributor distance from edge, m",
    )
    parser.add_argument(
        "--inlet-diameter", type=float, default=None, help="CF042 port diameter, m"
    )
    parser.add_argument(
        "--chamber-height",
        type=float,
        default=None,
        help="CF042 total compartment height, m",
    )
    parser.add_argument("--draw-mass-fraction", type=float, default=0.065)
    parser.add_argument("--feed-mass-fraction", type=float, default=0.00065)
    parser.add_argument(
        "--resistance",
        type=float,
        default=150666,
        help="2016 membrane support resistance, s/m",
    )
    parser.add_argument(
        "--results", type=Path, help="Export stem; default: CASE/result"
    )
    parser.add_argument(
        "--initialize-from",
        type=Path,
        help="Initialize CF042 from another mesh at the same geometry and operating point",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--setup-only", action="store_true", help="Prepare and mesh without solving"
    )
    action.add_argument(
        "--generate-only",
        action="store_true",
        help="Write case dictionaries without OpenFOAM",
    )
    action.add_argument(
        "--collect-only",
        action="store_true",
        help="Export results from an existing solved case",
    )
    parser.add_argument(
        "--momentum-relaxation",
        type=float,
        help="Transient momentum equation relaxation in (0, 1] (default: 1)",
    )
    args = parser.parse_args()
    if args.initialize_from and (
        args.example != "cf042" or args.generate_only or args.collect_only
    ):
        parser.error("--initialize-from requires CF042 mesh preparation or execution")
    if args.momentum_relaxation is not None and (
        not math.isfinite(args.momentum_relaxation)
        or not 0 < args.momentum_relaxation <= 1
    ):
        parser.error("--momentum-relaxation must be in (0, 1]")
    if args.transient is not None:
        if (
            args.end is not None
            or args.time_step is None
            or not all(
                math.isfinite(v) and v > 0 for v in (args.transient, args.time_step)
            )
            or args.time_step > args.transient
        ):
            parser.error(
                "Use --transient with a positive --time-step no larger than the duration, without --end"
            )
        steps = args.transient / args.time_step
        if not math.isfinite(steps) or not math.isclose(
            steps, round(steps), rel_tol=1e-10
        ):
            parser.error("Transient duration must be a multiple of the fixed time step")
        if args.outer_correctors is not None and args.outer_correctors < 1:
            parser.error("--outer-correctors must be positive")
        if (
            args.non_orthogonal_correctors is not None
            and args.non_orthogonal_correctors < 0
        ):
            parser.error("--non-orthogonal-correctors must be nonnegative")
    elif any(
        value is not None
        for value in (
            args.time_step,
            args.outer_correctors,
            args.non_orthogonal_correctors,
            args.momentum_relaxation,
        )
    ):
        parser.error("Time-step, corrector, and relaxation options require --transient")
    if args.ranks < 1 or (args.end is not None and args.end < 1):
        parser.error("--ranks and --end must be positive")
    if args.example != "channel" and (
        args.mode != "FO"
        or (args.generate_only and args.example not in ("cf042", "3dblock"))
    ):
        parser.error(
            "Chambers support FO; --generate-only supports channel, cf042 and 3dblock"
        )
    if args.example != "3dblock" and (args.spacers or args.spacer_radius is not None):
        parser.error("Spacer options apply to 3dblock")
    geometry_defaults = dict(
        length=0.085,
        width=0.039,
        height=0.00225,
        spacing=0.0008,
        layers=24,
        offset=0.0045,
        inlet_diameter=0.0055,
        chamber_height=0.01925,
    )
    if args.example == "3dblock" and any(
        getattr(args, k) is not None for k in geometry_defaults
    ):
        parser.error(
            "3dBlock has fixed 80 × 40 × 2 mm channels; dimension options apply to CF042"
        )
    for key, value in geometry_defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    if args.refinement != 1 and args.example not in (
        "3dblock",
        "chamber-a",
        "chamber-b",
    ):
        parser.error("--refinement applies to 2012 chambers and 3dBlock")
    case = args.case.resolve()
    output = args.results or case / "result"
    if args.collect_only:
        collect(case, output)
        return
    if case.exists():
        parser.error("Choose a new case directory")
    if not args.generate_only and not os.environ.get("WM_PROJECT_DIR"):
        parser.error("Source the OpenFOAM environment first")
    if args.ranks > 1 and os.environ.get("WM_COMPILER") == "Nvidia-gpu":
        parser.error("GPU execution supports one rank per case")

    def configure_time():
        if args.transient is None:
            return
        control = case / "system/controlDict"
        text = control.read_text()
        for key, value in dict(
            application="pisoSaltTransport",
            startFrom="startTime",
            startTime=0,
            endTime=args.transient,
            deltaT=args.time_step,
            writeControl="timeStep",
            writeInterval=max(1, math.ceil(args.transient / args.time_step / 20)),
        ).items():
            text = re.sub(rf"\b{key}\s+[^;]+;", f"{key} {value};", text)
        control.write_text(text)
        solution = case / "system/fvSolution"
        settings = re.sub(
            r"nOuterCorrectors\s+\d+;",
            f"nOuterCorrectors {args.outer_correctors or 4};",
            solution.read_text(),
        )
        nonorth = (
            2 if args.non_orthogonal_correctors is None else args.non_orthogonal_correctors
        )
        settings = re.sub(
            r"(PISO\s*\{[^}]*nNonOrthogonalCorrectors\s+)\d+",
            rf"\g<1>{nonorth}",
            settings,
        )
        # Clear SIMPLE relaxation; retain only explicitly requested momentum damping.
        momentum_relaxation = (
            1 if args.momentum_relaxation is None else args.momentum_relaxation
        )
        settings = re.sub(
            r"relaxationFactors\s*\{(?:[^{}]|\{[^{}]*\})*\}",
            (
                "relaxationFactors {}"
                if momentum_relaxation == 1
                else f"relaxationFactors {{ equations {{ U {momentum_relaxation}; }} }}"
            ),
            settings,
        )
        solution.write_text(settings)
        metadata = case / "case.json"
        parameters = json.loads(metadata.read_text())
        parameters.update(
            run_type="transient",
            end_time_s=args.transient,
            time_step_s=args.time_step,
            n_outer_correctors=args.outer_correctors or 4,
            n_non_orthogonal_correctors=nonorth,
            momentum_relaxation=momentum_relaxation,
        )
        metadata.write_text(json.dumps(parameters, indent=2) + "\n")

    def tool(name, *options):
        with (case / ("log." + name)).open("w") as log:
            subprocess.run(
                [name, "-case", str(case), *options],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    if args.example == "cf042":
        from membranefoam.cf042 import generate as cf042

        cf042(
            case,
            flow=args.flow,
            length=args.length,
            width=args.width,
            height=args.height,
            spacing=args.spacing,
            layers=args.layers,
            offset=args.offset,
            diameter=args.inlet_diameter,
            chamber_height=args.chamber_height,
            draw=args.draw_mass_fraction,
            feed=args.feed_mass_fraction,
            resistance=args.resistance,
            end=args.end or 60000,
        )
        if args.generate_only:
            configure_time()
            return
        if args.initialize_from:
            source = args.initialize_from.resolve()
            target_parameters = json.loads((case / "case.json").read_text())
            source_parameters = json.loads((source / "case.json").read_text())
            for key in (
                "geometry",
                "flow_ml_min",
                "length_m",
                "width_m",
                "height_m",
                "inlet_offset_m",
                "inlet_diameter_m",
                "chamber_height_m",
                "draw_mass_fraction",
                "feed_mass_fraction",
                "support_resistance_s_m",
            ):
                if source_parameters.get(key) != target_parameters[key]:
                    parser.error(f"Initialization source differs in {key}")
        tool("blockMesh")
        tool("checkMesh")
        if "Mesh OK" not in (case / "log.checkMesh").read_text():
            raise RuntimeError("Mesh check failed")
        tool("setFields")
        if args.initialize_from:
            options = [str(source), "-consistent", "-sourceTime", "latestTime"]
            if (source / "processor0").is_dir():
                options.append("-parallelSource")
            tool("mapFields", *options)
            # Coincident membrane faces are ambiguous to geometric interpolation.
            # Reset boundary guesses; the coupled conditions recompute both sides.
            with (case / "log.initializeBoundary").open("w") as log:
                for field, value in (
                    (
                        "m_A",
                        f"uniform {(args.draw_mass_fraction + args.feed_mass_fraction) / 2}",
                    ),
                    ("U", "uniform (0 0 0)"),
                ):
                    subprocess.run(
                        [
                            "foamDictionary",
                            str(case / "0" / field),
                            "-entry",
                            "boundaryField.membrane.value",
                            "-set",
                            value,
                        ],
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
        else:
            initialize_flow(case)
    elif args.example == "3dblock":
        from membranefoam import block

        build = block.generate if args.generate_only else block.prepare
        build(
            case,
            inlets=args.inlets,
            angle=args.angle,
            flow=args.flow,
            refinement=args.refinement,
            spacers=args.spacers,
            spacer_radius=args.spacer_radius,
            draw=args.draw_mass_fraction,
            feed=args.feed_mass_fraction,
            resistance=args.resistance,
            end=args.end or 60000,
        )
        if args.generate_only:
            configure_time()
            return
    elif args.example == "channel":
        generate(case, nx=args.nx, nz=args.nz, mode=args.mode, end=args.end or 12000)
        if args.generate_only:
            configure_time()
            return
        tool("blockMesh")
        tool("checkMesh")
        if "Mesh OK" not in (case / "log.checkMesh").read_text():
            raise RuntimeError("Mesh check failed")
    else:
        prepare(
            case, args.example[-1], end=args.end or 24000, refinement=args.refinement
        )
    configure_time()
    if args.setup_only:
        return
    solver = (
        "pisoSaltTransport" if args.transient is not None else "simpleSaltTransport"
    )
    command = [solver, "-case", str(case)]
    if args.ranks > 1:
        (case / "system/decomposeParDict").write_text(header("decomposeParDict") + f"""
numberOfSubdomains {args.ranks};
method simple;
simpleCoeffs {{ n ({args.ranks} 1 1); delta 0.001; }}
constraints {{ membranePairs {{ type preserveBaffles; }} }}
""")
        tool("decomposePar")
        command = ["mpirun", "-np", str(args.ranks), *command, "-parallel"]
    with (
        (case / "log.solver").open("w") as log,
        (case / "log.time").open("w") as timing,
    ):
        subprocess.run(
            ["/usr/bin/time", "-p", *command], stdout=log, stderr=timing, check=True
        )
    collect(case, output)
    print(f"Results: {output.with_name(output.name + '.json')}")


if __name__ == "__main__":
    main()
