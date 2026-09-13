"""Wall-contact spacer geometry derived from MembraneSimKit's 3dBlock.

Source revision: 2c4344144dadc3f2b486432724468a91652bbe3a.
Triangle and truncated-cylinder boundaries retain the source construction.
The cylinder uses corrected equal-area sizing; conformal transitions match
adjacent channel partitions without projecting spacer/linker interfaces.
"""

import math


def _wall_profile(shape, placement, fraction=0.2):
    h = 0.001
    H = 2 * h
    area = fraction * H * H
    profiles = []

    def add(q, n, g, arcs=()):
        profiles.append(dict(quad=q, counts=n, grading=g, arcs=arcs))

    if shape == "triangle":
        side = math.sqrt(4 * area / math.sqrt(3))
        height = math.sqrt(3) * side / 2
        left = -h + 0.1 * H
        mid = left + side / 2
        right = left + side
        add([(-h, 0), (left, 0), (mid, height), (-h, height)], (5, 10), (1, 5))
        add([(right, 0), (h, 0), (h, height), (mid, height)], (8, 10), (1, 5))
        add([(-h, height), (mid, height), (mid, H), (-h, H)], (5, 10), (1, 1))
        add([(mid, height), (h, height), (h, H), (mid, H)], (8, 10), (1, 1))
        footprint = side
    elif shape == "cylinder":
        radius = math.sqrt(area / (0.75 * math.pi + 0.5))
        c = radius / math.sqrt(2)
        add(
            [(-h, 0), (-c, 0), (-c, 2 * c), (-h, 2 * c)],
            (6, 10),
            (1, 5),
            [(1, 2, (-radius, c))],
        )
        add(
            [(c, 0), (h, 0), (h, 2 * c), (c, 2 * c)],
            (6, 10),
            (1, 5),
            [(3, 0, (radius, c))],
        )
        add([(-h, 2 * c), (-c, 2 * c), (-c, H), (-h, H)], (6, 10), (1, 1))
        add(
            [(-c, 2 * c), (c, 2 * c), (c, H), (-c, H)],
            (6, 10),
            (1, 1),
            [(0, 1, (0, radius + c))],
        )
        add([(c, 2 * c), (h, 2 * c), (h, H), (c, H)], (6, 10), (1, 1))
        footprint = 2 * c
    else:
        raise ValueError(shape)
    if placement == "up":
        for i, p in enumerate(profiles):
            lower = i >= 2
            p["quad"] = [(x, H - z) for x, z in reversed(p["quad"])]
            p["arcs"] = [(3 - a, 3 - b, (x, H - z)) for a, b, (x, z) in p["arcs"]]
            p["counts"] = (
                p["counts"][0],
                (15 if lower else 6) if shape == "cylinder" else 10,
            )
            p["grading"] = (1, 5 if lower else 1)
    elif placement != "down":
        raise ValueError(placement)
    return profiles, footprint


def wall_blocks(shape, placement, fraction, side_z, vertical):
    profiles, footprint = _wall_profile(shape, placement, fraction)
    h = 0.001
    H = 0.002
    out = []
    circle_radius = (
        math.sqrt(fraction * H * H / (0.75 * math.pi + 0.5))
        if shape == "cylinder"
        else None
    )
    center_z = (
        (
            circle_radius / math.sqrt(2)
            if placement == "down"
            else H - circle_radius / math.sqrt(2)
        )
        if circle_radius
        else None
    )
    for p in profiles:
        q = list(p["quad"])
        for i, (x, z) in enumerate(q):
            if abs(abs(x) - h) < 1e-12 and 1e-12 < z < H - 1e-12:
                q[i] = (x, h)
        outer = [z for x, z in q if abs(abs(x) - h) < 1e-12]
        lower = min(outer) == 0 if outer else placement == "up"
        cuts = (
            sorted(set([0.0, 1.0] + [z / h for z in side_z if 0 < z < h]))
            if lower
            else [0.0, 1.0]
        )
        arcs = {tuple(sorted((a, b))): (a, b, mid) for a, b, mid in p["arcs"]}

        def interpolate(a, b, t):
            if tuple(sorted((a, b))) in arcs:
                # A vertical curved edge is part of the circular boundary.
                theta = math.atan2(q[a][1] - center_z, q[a][0])
                end = math.atan2(q[b][1] - center_z, q[b][0])
                delta = (end - theta + math.pi) % (2 * math.pi) - math.pi
                angle = theta + t * delta
                return circle_radius * math.cos(
                    angle
                ), center_z + circle_radius * math.sin(angle)
            return tuple(q[a][j] + t * (q[b][j] - q[a][j]) for j in (0, 1))

        for lo, hi in zip(cuts, cuts[1:]):
            quad = [
                interpolate(0, 3, lo),
                interpolate(1, 2, lo),
                interpolate(1, 2, hi),
                interpolate(0, 3, hi),
            ]
            edges = []
            for a, b, mid in p["arcs"]:
                key = set((a, b))
                if key == {0, 3}:
                    edges.append((0, 3, interpolate(0, 3, (lo + hi) / 2)))
                elif key == {1, 2}:
                    edges.append((1, 2, interpolate(1, 2, (lo + hi) / 2)))
                elif key == {0, 1} and lo == 0:
                    edges.append((0, 1, mid))
                elif key == {2, 3} and hi == 1:
                    edges.append((2, 3, mid))
            nt, gt = vertical(lo * h, hi * h) if lower else (5, 1)
            out.append(
                dict(
                    quad=quad, counts=(p["counts"][0], nt), grading=(1, gt), arcs=edges
                )
            )
    return out, footprint
