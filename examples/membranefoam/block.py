"""3dBlock chamber with angled rectangular ports, based on the archived quarter mesh."""

import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from .channel import generate as channel, header
from .chambers import configure_2016, initialize_flow


def generate(
    path,
    inlets=3,
    angle=45,
    flow=50,
    refinement=1,
    draw=0.065,
    feed=0.00065,
    resistance=150666,
    end=60000,
    spacers=0,
    spacer_radius=None,
):
    """Write one quarter; prepare() mirrors its length and membrane plane."""
    if not isinstance(inlets, int) or inlets < 1 or inlets > 19 or inlets % 2 != 1:
        raise ValueError("Use an odd inlet count between 1 and 19")
    if not (0 <= angle <= 90 and refinement >= 1 and math.isfinite(refinement)):
        raise ValueError("Angle must be 0–90 degrees and refinement finite and >= 1")
    if not (
        all(math.isfinite(v) for v in (feed, draw, flow, resistance))
        and 0 <= feed < draw <= 0.09
        and flow > 0
        and resistance >= 0
    ):
        raise ValueError("Invalid operating conditions")
    if not isinstance(spacers, int) or spacers < 0 or spacers > 18 or spacers % 2:
        raise ValueError("Use an even spacer count between 0 and 18")
    if spacers and (spacer_radius is None or not 0 < spacer_radius < 0.001):
        raise ValueError("Specify a spacer radius strictly between 0 and 1 mm")
    if not spacers and spacer_radius is not None:
        raise ValueError("A spacer radius requires a nonzero spacer count")
    path = Path(path)
    channel(path, end=end)
    configure_2016(path, flow, draw, feed, resistance, end)
    # Tighter pressure tolerances can stall GAMG on the graded inlet mesh.
    solution = path / "system/fvSolution"
    solution.write_text(
        solution.read_text().replace("tolerance 1e-12;", "tolerance 1e-11;")
    )
    count = (inlets + 1) // 2
    theta = math.radians(angle)
    direction = (-math.cos(theta), math.sin(theta))
    tangent = (math.sin(theta), math.cos(theta))
    pivot = (0.0015, 0.0005)
    radius = 0.002

    def inlet_point(offset, outer=False):
        origin = tuple(pivot[i] + offset * tangent[i] for i in (0, 1))
        q = (origin[0] - radius, origin[1])
        dot = sum(q[i] * direction[i] for i in (0, 1))
        distance = (
            0.0056
            if outer
            else -dot + math.sqrt(dot * dot + radius * radius - sum(v * v for v in q))
        )
        return tuple(round(origin[i] + distance * direction[i], 14) for i in (0, 1))

    # Put block boundaries at both inlet lips. Cutting an opening through a
    # nearby grid line otherwise leaves thin, highly skewed wall cells.
    lips = [inlet_point(v) for v in (-0.0005, 0.0005)]
    angles = sorted(
        set(
            round(t, 12)
            for t in [0, math.pi / 4, math.pi / 2]
            + [math.atan2(z, radius - x) for x, z in lips]
        )
    )

    def ring_point(t, outer):
        r = radius if outer else radius / 2 / max(math.cos(t), math.sin(t))
        return (round(radius - r * math.cos(t), 14), round(r * math.sin(t), 14))

    inner = [ring_point(t, False) for t in angles]
    outer = [ring_point(t, True) for t in angles]
    xs = sorted(set([radius / 2, radius] + [x for x, z in inner]))
    zs = sorted(set([0.0, radius / 2] + [z for x, z in inner]))
    spans = [
        (
            0 if j == 0 else j * 0.04 / (inlets + 1),
            0.0005 if j == 0 else j * 0.04 / (inlets + 1) + 0.001,
        )
        for j in range(count)
    ]
    ys = sorted(set([0.0, 0.02] + [y for span in spans for y in span]))
    vertices, blocks, arcs, ports = [], [], [], []
    cache = {}
    junctions = set()

    def vertex(x, y, z):
        key = tuple(round(v, 14) for v in (x, y, z))
        if key not in cache:
            cache[key] = len(vertices)
            vertices.append(key)
        return cache[key]

    def horizontal(lo, hi):
        return max(1, math.ceil((hi - lo) / 0.000125 - 1e-8)), 1

    def vertical(lo, hi):
        ratio = (1 + 99 * hi / 0.001) / (1 + 99 * lo / 0.001)
        return max(1, math.ceil(15 * math.log(ratio) / math.log(100) - 1e-8)), ratio

    def extrude(quad, ylo, yhi, n, g):
        low = [vertex(x, ylo, z) for x, z in quad]
        high = [vertex(x, yhi, z) for x, z in quad]
        v = [low[0], low[1], high[1], high[0], low[3], low[2], high[2], high[3]]
        blocks.append(
            (v, (n[0], max(1, round((yhi - ylo) / 0.0001)), n[1]), (g[0], 1, g[1]))
        )
        return low, high

    for ylo, yhi in zip(ys, ys[1:]):
        for xlo, xhi in zip(xs, xs[1:]):
            for zlo, zhi in zip(zs, zs[1:]):
                nx, gx = horizontal(xlo, xhi)
                nz, gz = vertical(zlo, zhi)
                extrude(
                    [(xlo, zlo), (xhi, zlo), (xhi, zhi), (xlo, zhi)],
                    ylo,
                    yhi,
                    (nx, nz),
                    (gx, gz),
                )
        for i, (tlo, thi) in enumerate(zip(angles, angles[1:])):
            if thi <= math.pi / 4 + 1e-10:
                nt, gt = vertical(inner[i][1], inner[i + 1][1])
            else:
                nt, gt = horizontal(inner[i][0], inner[i + 1][0])
            low, high = extrude(
                [outer[i], inner[i], inner[i + 1], outer[i + 1]],
                ylo,
                yhi,
                (5, nt),
                (1, gt),
            )
            junctions.add(tuple(sorted([low[0], low[3], high[3], high[0]])))
            middle = ring_point((tlo + thi) / 2, True)
            for ids, y in [(low, ylo), (high, yhi)]:
                arcs.append((ids[0], ids[3], (middle[0], y, middle[1])))

        def plain_channel(xlo, xhi, uniform=False):
            for zlo, zhi in zip(zs + [radius], (zs + [radius])[1:]):
                nz, gz = vertical(zlo, zhi) if zhi <= 0.001 else (5, 1)
                extrude(
                    [(xlo, zlo), (xhi, zlo), (xhi, zhi), (xlo, zhi)],
                    ylo,
                    yhi,
                    (max(1, math.ceil((xhi - xlo) / 0.00025)) if uniform else 50, nz),
                    (1 if uniform else 8, gz),
                )

        def cylinder(centre):
            # O-grid around a transverse cylinder. The box sides retain the
            # neighbouring channel's vertical partitions and cell counts.
            h = radius / 2
            side_z = zs + [radius]
            perimeter = (
                [(h, z - h) for z in side_z if z >= h]
                + [(0, h), (-h, h)]
                + [(-h, z - h) for z in reversed(side_z) if z <= h]
                + [(0, -h), (h, -h)]
                + [(h, z - h) for z in side_z if 0 < z < h]
            )
            for i, a in enumerate(perimeter):
                b = perimeter[(i + 1) % len(perimeter)]
                ta, tb = math.atan2(a[1], a[0]), math.atan2(b[1], b[0])
                while tb <= ta:
                    tb += 2 * math.pi
                inner_a = (spacer_radius * math.cos(ta), spacer_radius * math.sin(ta))
                inner_b = (spacer_radius * math.cos(tb), spacer_radius * math.sin(tb))
                middle_a = tuple((a[j] + inner_a[j]) / 2 for j in (0, 1))
                middle_b = tuple((b[j] + inner_b[j]) / 2 for j in (0, 1))
                if abs(a[0] - b[0]) < 1e-12:
                    zlo, zhi = sorted((a[1] + h, b[1] + h))
                    nt, gt = vertical(zlo, zhi) if zhi <= h + 1e-12 else (5, 1)
                    if b[1] < a[1]:
                        gt = 1 / gt
                else:
                    nt, gt = 5, 1
                inner_mid = (
                    spacer_radius * math.cos((ta + tb) / 2),
                    spacer_radius * math.sin((ta + tb) / 2),
                )
                middle_mid = tuple(
                    (inner_mid[j] + (a[j] + b[j]) / 2) / 2 for j in (0, 1)
                )
                for q, grading, arc_mid in [
                    ([inner_a, middle_a, middle_b, inner_b], 32, inner_mid),
                    ([middle_a, a, b, middle_b], 1 / 32, middle_mid),
                ]:
                    low, high = extrude(
                        [(centre + x, h + z) for x, z in q],
                        ylo,
                        yhi,
                        (6, nt),
                        (grading, gt),
                    )
                    for ids, y in [(low, ylo), (high, yhi)]:
                        arcs.append(
                            (ids[0], ids[3], (centre + arc_mid[0], y, h + arc_mid[1]))
                        )

        start = radius
        for j in range(spacers // 2):
            centre = 0.08 * (j + 1) / (spacers + 1)
            plain_channel(start, centre - radius / 2, uniform=True)
            cylinder(centre)
            start = centre + radius / 2
        plain_channel(start, 0.04, uniform=bool(spacers))

    for lo, hi in spans:
        v = list(range(len(vertices), len(vertices) + 8))
        ports.append(v)
        for end in (False, True):
            for offset, y in [(-0.0005, lo), (0.0005, lo), (0.0005, hi), (-0.0005, hi)]:
                x, z = inlet_point(offset, end)
                vertices.append((x, y, z))
        middle = ring_point(sum(math.atan2(z, radius - x) for x, z in lips) / 2, True)
        arcs.extend(
            [
                (v[0], v[1], (middle[0], lo, middle[1])),
                (v[2], v[3], (middle[0], hi, middle[1])),
            ]
        )
        blocks.append((v, (9, max(1, round((hi - lo) / 0.0001)), 20), (1, 1, 5)))
    # Adjacent transverse blocks reference the same arc.
    arcs = list({tuple(sorted((a, b))): (a, b, p) for a, b, p in arcs}.values())
    faces = defaultdict(list)
    for v, n, g in blocks:
        for f in [
            (v[0], v[3], v[2], v[1]),
            (v[4], v[5], v[6], v[7]),
            (v[0], v[1], v[5], v[4]),
            (v[3], v[7], v[6], v[2]),
            (v[0], v[4], v[7], v[3]),
            (v[1], v[2], v[6], v[5]),
        ]:
            faces[tuple(sorted(f))].append(f)
    port_bases = {tuple(sorted(v[:4])) for v in ports}
    port_ends = {tuple(sorted(v[4:])) for v in ports}
    patches = defaultdict(list)
    for key, ff in faces.items():
        if len(ff) != 1:
            continue
        f = ff[0]
        points = [vertices[i] for i in f]
        if key in port_bases:
            name = "portBase"
        elif key in port_ends:
            name = "ports"
        elif key in junctions:
            name = "inletWall"
        elif all(abs(p[1]) < 1e-12 for p in points):
            name = "sides"
        elif all(abs(p[2]) < 1e-12 for p in points):
            name = "membrane"
        else:
            name = "walls"
        patches[name].append(f)

    def vector(v):
        return "(" + " ".join(map(str, v)) + ")"

    text = (
        header("blockMeshDict")
        + "scale 1;\nmergeType topology;\nvertices (\n"
        + "\n".join(map(vector, vertices))
        + "\n);\nblocks (\n"
    )
    for v, n, g in blocks:
        text += (
            "hex "
            + vector(v)
            + " "
            + vector([math.ceil(i * refinement) for i in n])
            + " simpleGrading "
            + vector(g)
            + "\n"
        )
    text += (
        ");\nedges (\n"
        + "\n".join(f"arc {a} {b} {vector(p)}" for a, b, p in arcs)
        + "\n);\nboundary (\n"
    )
    for name, ff in patches.items():
        kind = (
            "symmetryPlane"
            if name == "sides"
            else "patch" if name == "ports" else "wall"
        )
        text += (
            name
            + " { type "
            + kind
            + "; faces ("
            + " ".join(map(vector, ff))
            + "); }\n"
        )
    text += ");\nmergePatchPairs ((inletWall portBase));\n"
    (path / "system/blockMeshDict").write_text(text)
    (path / "case.json").write_text(
        json.dumps(
            dict(
                geometry="3dBlock reconstruction",
                paper="10.1016/j.seppur.2015.12.017",
                length_m=0.08,
                width_m=0.04,
                height_m=0.002,
                inlet_count=inlets,
                inlet_angle_degrees=angle,
                inlet_width_m=0.001,
                flow_ml_min=flow,
                draw_mass_fraction=draw,
                feed_mass_fraction=feed,
                support_resistance_s_m=resistance,
                refinement=refinement,
                half_width=True,
                spacer_count=spacers,
                spacer_shape="cylinder" if spacers else None,
                spacer_placement="middle" if spacers else None,
                spacer_radius_m=spacer_radius,
            ),
            indent=2,
        )
        + "\n"
    )


def prepare(path, **parameters):
    path = Path(path).resolve()
    generate(path, **parameters)
    initial = path / "initialFields"
    (path / "0").rename(initial)

    def tool(name, *args):
        with (path / ("log." + name)).open("a") as log:
            subprocess.run(
                [name, "-case", str(path), *args],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    tool("blockMesh")
    for point, normal in [(".04 0 0", "1 0 0"), ("0 0 0", "0 0 1")]:
        (path / "system/mirrorMeshDict").write_text(header("mirrorMeshDict") + f"""
planeType pointAndNormal;
pointAndNormalDict {{ basePoint ({point}); normalVector ({normal}); }}
planeTolerance 1e-8;
""")
        tool("mirrorMesh", "-overwrite")
    actions = [
        "{ name membraneFaces; type faceSet; action new; source boxToFace; box (0 -1 -1e-10) (.08 1 1e-10); }",
        "{ name membraneZone; type faceZoneSet; action new; source setToFaceZone; faceSet membraneFaces; }",
        "{ name portFloorFaces; type faceSet; action new; source boxToFace; box (-1 -1 -1e-10) (1 1 1e-10); }",
        "{ name portFloorFaces; type faceSet; action subtract; source boxToFace; box (0 -1 -1e-10) (.08 1 1e-10); }",
        "{ name portFloorZone; type faceZoneSet; action new; source setToFaceZone; faceSet portFloorFaces; }",
    ]
    patch_defs = []
    for name, xlo, xhi, zlo, zhi in [
        ("feedIn", -1, 0.04, -1, 0),
        ("feedOut", 0.04, 1, -1, 0),
        ("drawIn", 0.04, 1, 0, 1),
        ("drawOut", -1, 0.04, 0, 1),
    ]:
        actions += [
            f"{{ name {name}; type faceSet; action new; source patchToFace; patches (ports); }}",
            f"{{ name {name}; type faceSet; action subset; source boxToFace; box ({xlo} -1 {zlo}) ({xhi} 1 {zhi}); }}",
        ]
        patch_defs.append(
            f"{{ name {name}; patchInfo {{ type patch; }} constructFrom set; set {name}; }}"
        )
    patch_defs.append(
        "{ name walls; patchInfo { type wall; } constructFrom patches; patches (walls inletWall portBase); }"
    )
    (path / "system/topoSetDict").write_text(
        header("topoSetDict") + "actions (" + "\n".join(actions) + ");\n"
    )
    tool("topoSet")
    (path / "system/createBafflesDict").write_text(header("createBafflesDict") + """
internalFacesOnly true;
baffles { membrane { type faceZone; zoneName membraneZone;
 patches { master { name membrane; type wall; } slave { name membrane; type wall; } }
}
portFloors { type faceZone; zoneName portFloorZone;
 patches { master { name walls; type wall; } slave { name walls; type wall; } }
} }
""")
    tool("createBaffles", "-overwrite")
    # Re-select port faces after baffle creation changes face numbering.
    tool("topoSet")
    (path / "system/createPatchDict").write_text(
        header("createPatchDict")
        + "pointSync false; patches ("
        + "\n".join(patch_defs)
        + ");\n"
    )
    tool("createPatch", "-overwrite")
    tool("checkMesh")
    if "Mesh OK" not in (path / "log.checkMesh").read_text():
        raise ValueError("Mesh check failed")
    (path / "0").mkdir(exist_ok=True)
    for field in initial.iterdir():
        field.rename(path / "0" / field.name)
    initial.rmdir()
    tool("setFields")
    initialize_flow(path)
