"""Structured reconstruction of the CF042 chamber used in the 2016 study.

Coordinates are SI. The membrane is z=0; y=0 is the half-width symmetry plane.
Rounded distributors connect the thin channels to circular, normal inlet ports.
"""

import json
import math
from collections import defaultdict
from pathlib import Path
from .channel import generate as channel, header
from .chambers import configure_2016


def generate(
    path,
    flow=50,
    length=0.085,
    width=0.039,
    height=0.00225,
    spacing=0.0008,
    layers=24,
    offset=0.0045,
    diameter=0.0055,
    chamber_height=0.01925,
    draw=0.065,
    feed=0.00065,
    resistance=150666.0,
    end=60000,
):
    if not all(
        math.isfinite(value)
        for value in (
            flow,
            length,
            width,
            height,
            spacing,
            layers,
            offset,
            diameter,
            chamber_height,
            draw,
            feed,
            resistance,
            end,
        )
    ):
        raise ValueError("All chamber parameters must be finite")
    r = diameter / 2
    half = width / 2
    span = half - r
    plenum_top = chamber_height - r
    if not (
        0 <= feed < draw <= 0.09
        and flow > 0
        and layers >= 4
        and int(layers) == layers
        and spacing > 0
        and resistance >= 0
    ):
        raise ValueError(
            "Require 0 <= feed < draw <= 0.09, positive flow/spacing, "
            "integer layers >= 4 and nonnegative resistance"
        )
    if not (
        0 < r
        and 0 < offset
        and 2 * (offset + diameter) < length
        and half > 4 * r
        and 0 < height < plenum_top
    ):
        raise ValueError("Chamber dimensions do not leave room for inlet distributors")
    path = Path(path)
    channel(path, end=end)
    xs = sorted(
        set(
            [
                0.0,
                r,
                offset,
                offset + r,
                offset + diameter,
                length - offset - diameter,
                length - offset - r,
                length - offset,
                length - r,
                length,
            ]
        )
    )
    ys = [0.0, r, span - r, span, half - r, half]
    if offset < r:
        ys += [span - r / 2, half - r / 2]
    ys = sorted(set(ys))
    zs = [0.0, height, plenum_top, chamber_height]
    centres = [offset + r, length - offset - r]
    vertices = []
    index = {}
    blocks = []
    faces = defaultdict(list)
    arcs = {}

    def xy(x, y, z):
        # Square-to-disk mapping keeps intermediate grid lines ordered along arcs.
        for c, sign in [(r, -1), (length - r, 1)]:
            u = sign * (x - c) / r
            v = (y - half + r) / r
            if 0 <= u <= 1 + 1e-12 and 0 < v <= 1 + 1e-12:
                return c + sign * r * u * math.sqrt(
                    max(0, 1 - v * v / 2)
                ), half - r + r * v * math.sqrt(max(0, 1 - u * u / 2))
        for c in centres:
            u = (x - c) / r
            v = (y - span + r) / r
            if abs(u) <= 1 + 1e-12 and 0 < v <= 1 + 1e-12:
                return c + r * u * math.sqrt(
                    max(0, 1 - v * v / 2)
                ), span - r + r * v * math.sqrt(max(0, 1 - u * u / 2))
            v = y / r
            if z >= plenum_top and abs(u) <= 1 + 1e-12 and 0 <= v <= 1 + 1e-12:
                return c + r * u * math.sqrt(max(0, 1 - v * v / 2)), r * v * math.sqrt(
                    max(0, 1 - u * u / 2)
                )
        return x, y

    def vertex(side, i, j, k):
        key = side, i, j, k
        if key not in index:
            x, y = xy(xs[i], ys[j], zs[k])
            index[key] = len(vertices)
            vertices.append((x, y, side * zs[k]))
        return index[key]

    def arc(a, b, cx, cy):
        if tuple(sorted((a, b))) in arcs:
            return
        pa, pb = vertices[a], vertices[b]
        va = (pa[0] - cx, pa[1] - cy)
        vb = (pb[0] - cx, pb[1] - cy)
        vv = (va[0] + vb[0], va[1] + vb[1])
        norm = math.hypot(*vv)
        arcs[tuple(sorted((a, b)))] = (
            cx + r * vv[0] / norm,
            cy + r * vv[1] / norm,
            pa[2],
        )

    for side in [-1, 1]:
        for i in range(len(xs) - 1):
            for j in range(len(ys) - 1):
                xm = (xs[i] + xs[i + 1]) / 2
                ym = (ys[j] + ys[j + 1]) / 2
                distributor = any(abs(xm - c) < r for c in centres) and ym < span
                pipe = distributor and ym < r
                for k in range(3):
                    if k == 1 and not distributor or k == 2 and not pipe:
                        continue
                    lo, hi = (k, k + 1) if side == 1 else (k + 1, k)
                    v = [
                        vertex(side, i, j, lo),
                        vertex(side, i + 1, j, lo),
                        vertex(side, i + 1, j + 1, lo),
                        vertex(side, i, j + 1, lo),
                        vertex(side, i, j, hi),
                        vertex(side, i + 1, j, hi),
                        vertex(side, i + 1, j + 1, hi),
                        vertex(side, i, j + 1, hi),
                    ]
                    nx = max(2, math.ceil((xs[i + 1] - xs[i]) / spacing))
                    ny = max(2, math.ceil((ys[j + 1] - ys[j]) / spacing))
                    nz = (
                        layers
                        if k == 0
                        else max(3, math.ceil((zs[k + 1] - zs[k]) / spacing))
                    )
                    grading = (100 if side == 1 else 0.01) if k == 0 else 1
                    blocks.append(
                        f"hex ({' '.join(map(str, v))}) ({nx} {ny} {nz}) simpleGrading (1 1 {grading})"
                    )
                    for f in [
                        (v[0], v[3], v[2], v[1]),
                        (v[4], v[5], v[6], v[7]),
                        (v[0], v[1], v[5], v[4]),
                        (v[3], v[7], v[6], v[2]),
                        (v[0], v[4], v[7], v[3]),
                        (v[1], v[2], v[6], v[5]),
                    ]:
                        faces[tuple(sorted(f))].append(f)
                    # Curved horizontal edges on each extrusion plane.
                    for kk in (k, k + 1):
                        for (ia, ja), (ib, jb) in [
                            ((i, j), (i + 1, j)),
                            ((i + 1, j), (i + 1, j + 1)),
                            ((i + 1, j + 1), (i, j + 1)),
                            ((i, j + 1), (i, j)),
                        ]:
                            a, b = vertex(side, ia, ja, kk), vertex(side, ib, jb, kk)
                            pa, pb = vertices[a], vertices[b]
                            circles = []
                            if k == 0:
                                circles += [(r, half - r), (length - r, half - r)]
                            circles += [(c, span - r) for c in centres]
                            if zs[kk] >= plenum_top:
                                circles += [(c, 0) for c in centres]
                            for cx, cy in circles:
                                if (
                                    abs(math.hypot(pa[0] - cx, pa[1] - cy) - r) < 1e-10
                                    and abs(math.hypot(pb[0] - cx, pb[1] - cy) - r)
                                    < 1e-10
                                ):
                                    arc(a, b, cx, cy)
    patches = defaultdict(list)
    for vv in faces.values():
        if len(vv) != 1:
            continue
        f = vv[0]
        pts = [vertices[v] for v in f]
        if all(abs(p[2]) < 1e-12 for p in pts):
            name = "membrane"
        elif all(abs(p[1]) < 1e-12 for p in pts):
            name = "sides"
        elif all(abs(abs(p[2]) - chamber_height) < 1e-12 for p in pts):
            feed_side = pts[0][2] < 0
            left = sum(p[0] for p in pts) / 4 < length / 2
            name = ("feed" if feed_side else "draw") + (
                "In" if left == feed_side else "Out"
            )
        else:
            name = "walls"
        patches[name].append("(" + " ".join(map(str, f)) + ")")
    text = header("blockMeshDict") + "scale 1; mergeType topology;\nvertices (\n"
    text += (
        "\n".join(f"({x:.14g} {y:.14g} {z:.14g})" for x, y, z in vertices)
        + "\n);\nblocks (\n"
        + "\n".join(blocks)
        + "\n);\nedges (\n"
    )
    text += (
        "\n".join(
            f"arc {a} {b} ({x:.14g} {y:.14g} {z:.14g})"
            for (a, b), (x, y, z) in arcs.items()
        )
        + "\n);\nboundary (\n"
    )
    for name, ff in patches.items():
        kind = (
            "symmetryPlane"
            if name == "sides"
            else "wall"
            if name in ("walls", "membrane")
            else "patch"
        )
        text += f"{name} {{ type {kind}; faces ( {' '.join(ff)} ); }}\n"
    (path / "system/blockMeshDict").write_text(text + ");\nmergePatchPairs ();\n")
    configure_2016(path, flow, draw, feed, resistance, end)
    (path / "case.json").write_text(
        json.dumps(
            dict(
                geometry="CF042 reconstruction",
                paper="10.1016/j.seppur.2015.12.017",
                flow_ml_min=flow,
                length_m=length,
                width_m=width,
                height_m=height,
                inlet_offset_m=offset,
                inlet_diameter_m=diameter,
                chamber_height_m=chamber_height,
                draw_mass_fraction=draw,
                feed_mass_fraction=feed,
                support_resistance_s_m=resistance,
                spacing_m=spacing,
                membrane_layers=layers,
                wall_normal_expansion=100,
                channel_corner_radius_m=r,
                distributor_half_span_m=span,
                plenum_top_m=plenum_top,
                half_width=True,
            ),
            indent=2,
        )
        + "\n"
    )
