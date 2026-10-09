"""Recovered CF042 block construction with current mesh and field preparation.

Geometry derives from Mathias Gruber's MembraneSimKit blockCF042 generator,
commit 2c4344144dadc3f2b486432724468a91652bbe3a. Mesh grading and topology are
retained; duplicate identical arcs are removed for OpenFOAM v2606.
"""

import math
from math import floor
from collections import Counter
from .channel import header


class QuarterMesh:
    """Build the original quarter mesh without the legacy Python-2 toolkit."""

    def __init__(self):
        self.points = []
        self.indices = {}
        self.blocks = []
        self.edges = {}
        self.patches = {}
        self.merge_pairs = []

    def patch(self, name, kind):
        self.patches.setdefault(name, dict(kind=kind, faces=[]))

    def merge(self, master, slave):
        self.merge_pairs.append((master, slave))

    @staticmethod
    def cos_degrees(angle):
        return math.cos(math.radians(angle))

    @staticmethod
    def sin_degrees(angle):
        return math.sin(math.radians(angle))

    def square(self, first, second, names, counts, grading):
        lo = [min(a, b) for a, b in zip(first, second)]
        hi = [max(a, b) for a, b in zip(first, second)]
        return self.custom(
            [
                [x, y, z]
                for z in (lo[2], hi[2])
                for x, y in (
                    (lo[0], lo[1]),
                    (hi[0], lo[1]),
                    (hi[0], hi[1]),
                    (lo[0], hi[1]),
                )
            ],
            names,
            counts,
            grading,
        )

    def custom(self, points, names, counts, grading):
        ids = []
        for point in points:
            key = tuple(point)
            if key not in self.indices:
                self.indices[key] = len(self.points)
                self.points.append(tuple(point))
            ids.append(self.indices[key])
        counts = tuple(max(1, round(n)) for n in counts)
        self.blocks.append((ids, counts, grading))
        for name, face in zip(
            names,
            [
                (0, 3, 2, 1),
                (0, 1, 5, 4),
                (0, 4, 7, 3),
                (6, 5, 1, 2),
                (3, 7, 6, 2),
                (4, 5, 6, 7),
            ],
        ):
            self.patch(name, "patch")
            self.patches[name]["faces"].append(tuple(ids[i] for i in face))
        return len(self.blocks) - 1, ids

    def arc(self, a, b, midpoint):
        key = tuple(sorted((a, b)))
        if key in self.edges and self.edges[key][2] != tuple(midpoint):
            raise ValueError(f"Conflicting curved edge {key}")
        self.edges[key] = a, b, tuple(midpoint)

    def dictionary(self):
        def values(row):
            return "(" + " ".join(str(v) for v in row) + ")"

        text = header("blockMeshDict") + "scale 1;\nvertices (\n"
        text += "\n".join(values(p) for p in self.points) + "\n);\nblocks (\n"
        text += (
            "\n".join(
                f"hex {values(v)} {values(n)} simpleGrading {values(g)}"
                for v, n, g in self.blocks
            )
            + "\n);\nedges (\n"
        )
        text += "\n".join(f"arc {a} {b} {values(p)}" for a, b, p in self.edges.values())
        text += "\n);\nboundary (\n"
        counts = Counter(
            tuple(sorted(f)) for p in self.patches.values() for f in p["faces"]
        )
        for name, patch in self.patches.items():
            faces = [f for f in patch["faces"] if counts[tuple(sorted(f))] == 1]
            if faces:
                text += f"{name} {{ type {patch['kind']}; faces (\n"
                text += "\n".join(values(f) for f in faces) + "\n); }\n"
        text += (
            ");\nmergePatchPairs ("
            + " ".join(values(p) for p in self.merge_pairs)
            + ");\n"
        )
        return text


class Geometry:

    def __init__(
        self,
        *,
        length,
        width,
        height,
        diameter,
        corner_radius,
        offset,
        distributor_width,
        inlet_height,
        refinement,
        membrane_expansion=5,
    ):
        inletDia = diameter
        cornerRad = corner_radius
        edgeDist = offset
        inletwidth = distributor_width
        inletheight = inlet_height
        fine = refinement
        mesh = QuarterMesh()
        mesh.patch("chamberWall", "wall")
        mesh.patch("inletWall", "wall")
        mesh.patch("tempWall1", "wall")
        mesh.patch("tempWall2", "wall")
        mesh.patch("inlet", "patch")
        mesh.patch("outlet", "patch")
        mesh.patch("membrane", "patch")
        mesh.patch("symm", "symmetryPlane")
        self.h1 = height
        self.h2 = height / 2
        self.w2 = width / 2
        self.iw2 = inletwidth / 2
        self.iw4 = inletwidth / 4
        self.l2 = length / 2
        self.ID = inletDia
        self.IR = inletDia / 2
        self.R = cornerRad
        self.R2 = self.R / 2
        self.cos45 = mesh.cos_degrees(45)
        self.cos22half = mesh.cos_degrees(22.5)
        self.sin22half = mesh.sin_degrees(22.5)
        self.fine = fine
        self.membraneScale = membrane_expansion
        self.symmScale = 1
        self.iScale = 0.5
        self.outToIn = (1 - self.iScale) * self.R
        self.bufferDist = 0.0005
        self.inletStart = edgeDist + self.bufferDist
        self.mergeTol = 1e-06
        self.xCells = int(length / 0.085 * 150)
        self.yCells = int((self.w2 - self.R - self.iw2) * 100 * 30)
        h1 = self.h1
        w2 = self.w2
        iw2 = self.iw2
        iw4 = self.iw4
        l2 = self.l2
        IR = self.IR
        R = self.R
        outToIn = self.outToIn
        cos45 = self.cos45
        cos22half = self.cos22half
        sin22half = self.sin22half
        membraneScale = self.membraneScale
        symmScale = self.symmScale
        inletStart = self.inletStart
        mergeTol = self.mergeTol
        xCells = self.xCells
        yCells = self.yCells
        xMid, yMid = (R - mesh.cos_degrees(45) * R, w2 - (R - mesh.cos_degrees(45) * R))
        mesh.square(
            [outToIn, w2 - R, 0],
            [R, w2 - outToIn, h1],
            [
                "membrane",
                "tempWall1",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 5, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        block1, points1 = mesh.custom(
            [
                [0, w2 - R, 0],
                [outToIn, w2 - R, 0],
                [outToIn, w2 - outToIn, 0],
                [xMid, yMid, 0],
                [0, w2 - R, h1],
                [outToIn, w2 - R, h1],
                [outToIn, w2 - outToIn, h1],
                [xMid, yMid, h1],
            ],
            [
                "membrane",
                "tempWall1",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 5, fine * 5, fine * 25],
            [1, symmScale, membraneScale],
        )
        mesh.arc(
            points1[0], points1[3], [R - cos22half * R, w2 - (R - sin22half * R), 0]
        )
        mesh.arc(
            points1[4], points1[7], [R - cos22half * R, w2 - (R - sin22half * R), h1]
        )
        block2, points2 = mesh.custom(
            [
                [outToIn, w2 - outToIn, 0],
                [R, w2 - outToIn, 0],
                [R, w2, 0],
                [xMid, yMid, 0],
                [outToIn, w2 - outToIn, h1],
                [R, w2 - outToIn, h1],
                [R, w2, h1],
                [xMid, yMid, h1],
            ],
            [
                "membrane",
                "tempWall1",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 5, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.arc(
            points2[2], points2[3], [R - sin22half * R, w2 - (R - cos22half * R), 0]
        )
        mesh.arc(
            points2[6], points2[7], [R - sin22half * R, w2 - (R - cos22half * R), h1]
        )
        mesh.square(
            [R, w2 - outToIn, 0],
            [l2, w2, h1],
            [
                "membrane",
                "tempWall1",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * xCells, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [R, w2 - R, 0],
            [l2, w2 - outToIn, h1],
            [
                "membrane",
                "tempWall1",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * xCells, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        self.inlet_region(mesh, 0, h1, 25, membraneScale)
        mesh.custom(
            [
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, 0],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, 0],
                [inletStart + R + cos45 * R, w2 - R - mergeTol, 0],
                [inletStart + R - cos45 * R, w2 - R - mergeTol, 0],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, h1],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, h1],
                [inletStart + R + cos45 * R, w2 - R - mergeTol, h1],
                [inletStart + R - cos45 * R, w2 - R - mergeTol, h1],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "tempWall2",
                "chamberWall",
            ],
            [fine * 12, fine * yCells, fine * 25],
            [1, symmScale, membraneScale],
        )
        mesh.merge("tempWall2", "tempWall1")
        fraction = edgeDist / (l2 - 2 * R)
        leftSideCells = int(floor(xCells * fraction))
        mesh.square(
            [0, 0, 0],
            [inletStart, cos45 * IR, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * leftSideCells, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [0, cos45 * IR, 0],
            [inletStart, iw4, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * leftSideCells, fine * 10, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [0, iw4, 0],
            [inletStart, iw2 - R, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * leftSideCells, fine * 10, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [0, iw2 - R + cos45 * R, 0],
            [inletStart + R - cos45 * R, w2 - R - mergeTol, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "tempWall2",
                "chamberWall",
            ],
            [fine * leftSideCells, fine * yCells, fine * 25],
            [1, 1, membraneScale],
        )
        inletLeftBlock, inletLeftPoints = mesh.custom(
            [
                [0, iw2 - R, 0],
                [inletStart, iw2 - R, 0],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, 0],
                [0, iw2 - R + cos45 * R, 0],
                [0, iw2 - R, h1],
                [inletStart, iw2 - R, h1],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, h1],
                [0, iw2 - R + cos45 * R, h1],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * leftSideCells, fine * 5, fine * 25],
            [1, symmScale, membraneScale],
        )
        rightSideCells = int(floor(xCells * (1 - fraction)))
        mesh.square(
            [inletStart + 2 * R, 0, 0],
            [l2, cos45 * IR, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * rightSideCells, fine * 5, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [inletStart + 2 * R, cos45 * IR, 0],
            [l2, iw4, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * rightSideCells, fine * 10, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [inletStart + 2 * R, iw4, 0],
            [l2, iw2 - R, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * rightSideCells, fine * 10, fine * 25],
            [1, 1, membraneScale],
        )
        mesh.square(
            [inletStart + R + cos45 * R, iw2 - R + cos45 * R, 0],
            [l2, w2 - R - mergeTol, h1],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "tempWall2",
                "chamberWall",
            ],
            [fine * rightSideCells, fine * yCells, fine * 25],
            [1, 1, membraneScale],
        )
        inletLeftBlock, inletLeftPoints = mesh.custom(
            [
                [inletStart + 2 * R, iw2 - R, 0],
                [l2, iw2 - R, 0],
                [l2, iw2 - R + cos45 * R, 0],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, 0],
                [inletStart + 2 * R, iw2 - R, h1],
                [l2, iw2 - R, h1],
                [l2, iw2 - R + cos45 * R, h1],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, h1],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * rightSideCells, fine * 5, fine * 25],
            [1, symmScale, membraneScale],
        )
        self.inlet_region(mesh, h1, inletheight, 25, 1)
        self.inlet_tube(mesh, inletheight, inletheight + h1, 6, 1)
        self.mesh = mesh

    def inlet_region(self, mesh, lowZ, highZ, cells, scaling):
        iw2 = self.iw2
        iw4 = self.iw4
        IR = self.IR
        R = self.R
        outToIn = self.outToIn
        fine = self.fine
        cos45 = self.cos45
        cos22half = self.cos22half
        sin22half = self.sin22half
        symmScale = self.symmScale
        iScale = self.iScale
        inletStart = self.inletStart
        self.inlet_tube(mesh, lowZ, highZ, cells, scaling)
        inletRightMostBlock, inletRightMostPoints = mesh.custom(
            [
                [inletStart + R + IR, 0, lowZ],
                [inletStart + 2 * R, 0, lowZ],
                [inletStart + 2 * R, cos45 * IR, lowZ],
                [inletStart + R + cos45 * IR, cos45 * IR, lowZ],
                [inletStart + R + IR, 0, highZ],
                [inletStart + 2 * R, 0, highZ],
                [inletStart + 2 * R, cos45 * IR, highZ],
                [inletStart + R + cos45 * IR, cos45 * IR, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        inletLeftMostBlock, inletLeftMostPoints = mesh.custom(
            [
                [inletStart, 0, lowZ],
                [inletStart + R - IR, 0, lowZ],
                [inletStart + R - IR * cos45, cos45 * IR, lowZ],
                [inletStart, cos45 * IR, lowZ],
                [inletStart, 0, highZ],
                [inletStart + R - IR, 0, highZ],
                [inletStart + R - IR * cos45, cos45 * IR, highZ],
                [inletStart, cos45 * IR, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.square(
            [inletStart, cos45 * IR, lowZ],
            [inletStart + R - IR * cos45, iw4, highZ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 10, fine * cells],
            [1, 1, scaling],
        )
        mesh.square(
            [inletStart + R - IR * cos45, cos45 * IR, lowZ],
            [inletStart + R + IR * cos45, iw4, highZ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 12, fine * 10, fine * cells],
            [1, 1, scaling],
        )
        mesh.square(
            [inletStart + 2 * R, cos45 * IR, lowZ],
            [inletStart + R + IR * cos45, iw4, highZ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 10, fine * cells],
            [1, 1, scaling],
        )
        mesh.custom(
            [
                [inletStart + R - cos45 * IR, iw4, lowZ],
                [inletStart + R + cos45 * IR, iw4, lowZ],
                [inletStart + R + iScale * R, iw2 - R, lowZ],
                [inletStart + R - iScale * R, iw2 - R, lowZ],
                [inletStart + R - cos45 * IR, iw4, highZ],
                [inletStart + R + cos45 * IR, iw4, highZ],
                [inletStart + R + iScale * R, iw2 - R, highZ],
                [inletStart + R - iScale * R, iw2 - R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 12, fine * 10, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.custom(
            [
                [inletStart, iw4, lowZ],
                [inletStart + R - cos45 * IR, iw4, lowZ],
                [inletStart + R - iScale * R, iw2 - R, lowZ],
                [inletStart, iw2 - R, lowZ],
                [inletStart, iw4, highZ],
                [inletStart + R - cos45 * IR, iw4, highZ],
                [inletStart + R - iScale * R, iw2 - R, highZ],
                [inletStart, iw2 - R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 10, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.custom(
            [
                [inletStart + R + cos45 * IR, iw4, lowZ],
                [inletStart + 2 * R, iw4, lowZ],
                [inletStart + 2 * R, iw2 - R, lowZ],
                [inletStart + R + iScale * R, iw2 - R, lowZ],
                [inletStart + R + cos45 * IR, iw4, highZ],
                [inletStart + 2 * R, iw4, highZ],
                [inletStart + 2 * R, iw2 - R, highZ],
                [inletStart + R + iScale * R, iw2 - R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 10, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.square(
            [inletStart + R - iScale * R, iw2 - R, lowZ],
            [inletStart + R + iScale * R, iw2 - R + iScale * R, highZ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 12, fine * 5, fine * cells],
            [1, 1, scaling],
        )
        inletLeftBlock, inletLeftPoints = mesh.custom(
            [
                [inletStart, iw2 - R, lowZ],
                [inletStart + outToIn, iw2 - R, lowZ],
                [inletStart + outToIn, iw2 - outToIn, lowZ],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, lowZ],
                [inletStart, iw2 - R, highZ],
                [inletStart + outToIn, iw2 - R, highZ],
                [inletStart + outToIn, iw2 - outToIn, highZ],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(
            inletLeftPoints[0],
            inletLeftPoints[3],
            [inletStart + R - cos22half * R, iw2 - R + sin22half * R, lowZ],
        )
        mesh.arc(
            inletLeftPoints[4],
            inletLeftPoints[7],
            [inletStart + R - cos22half * R, iw2 - R + sin22half * R, highZ],
        )
        inletRightBlock, inletRightPoints = mesh.custom(
            [
                [inletStart + 2 * R - outToIn, iw2 - R, lowZ],
                [inletStart + 2 * R, iw2 - R, lowZ],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, lowZ],
                [inletStart + 2 * R - outToIn, iw2 - R + iScale * R, lowZ],
                [inletStart + 2 * R - outToIn, iw2 - R, highZ],
                [inletStart + 2 * R, iw2 - R, highZ],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, highZ],
                [inletStart + 2 * R - outToIn, iw2 - R + iScale * R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 3, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(
            inletRightPoints[1],
            inletRightPoints[2],
            [inletStart + R + cos22half * R, iw2 - R + sin22half * R, lowZ],
        )
        mesh.arc(
            inletRightPoints[5],
            inletRightPoints[6],
            [inletStart + R + cos22half * R, iw2 - R + sin22half * R, highZ],
        )
        inletTopBlock, inletTopPoints = mesh.custom(
            [
                [inletStart + R - iScale * R, iw2 - outToIn, lowZ],
                [inletStart + R + iScale * R, iw2 - outToIn, lowZ],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, lowZ],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, lowZ],
                [inletStart + R - iScale * R, iw2 - outToIn, highZ],
                [inletStart + R + iScale * R, iw2 - outToIn, highZ],
                [inletStart + R + cos45 * R, iw2 - R + cos45 * R, highZ],
                [inletStart + R - cos45 * R, iw2 - R + cos45 * R, highZ],
            ],
            [
                "membrane",
                "symm",
                "chamberWall",
                "chamberWall",
                "chamberWall",
                "chamberWall",
            ],
            [fine * 12, fine * 3, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(inletTopPoints[2], inletTopPoints[3], [inletStart + R, iw2, lowZ])
        mesh.arc(inletTopPoints[6], inletTopPoints[7], [inletStart + R, iw2, highZ])

    def inlet_tube(self, mesh, lowZ, highZ, cells, scaling):
        IR = self.IR
        R = self.R
        fine = self.fine
        cos45 = self.cos45
        cos22half = self.cos22half
        sin22half = self.sin22half
        symmScale = self.symmScale
        iScale = self.iScale
        inletStart = self.inletStart
        inletCenterBlock, inletCenterPoints = mesh.square(
            [inletStart + R - iScale * IR, 0, lowZ],
            [inletStart + R + iScale * IR, iScale * IR, highZ],
            ["membrane", "symm", "chamberWall", "chamberWall", "chamberWall", "inlet"],
            [fine * 12, fine * 5, fine * cells],
            [1, 1, scaling],
        )
        inletLeftBlock, inletLeftPoints = mesh.custom(
            [
                [inletStart + R - IR, 0, lowZ],
                [inletStart + R - iScale * IR, 0, lowZ],
                [inletStart + R - iScale * IR, iScale * IR, lowZ],
                [inletStart + R - cos45 * IR, cos45 * IR, lowZ],
                [inletStart + R - IR, 0, highZ],
                [inletStart + R - iScale * IR, 0, highZ],
                [inletStart + R - iScale * IR, iScale * IR, highZ],
                [inletStart + R - cos45 * IR, cos45 * IR, highZ],
            ],
            ["membrane", "symm", "chamberWall", "chamberWall", "chamberWall", "inlet"],
            [fine * 5, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(
            inletLeftPoints[0],
            inletLeftPoints[3],
            [inletStart + R - cos22half * IR, sin22half * IR, lowZ],
        )
        mesh.arc(
            inletLeftPoints[4],
            inletLeftPoints[7],
            [inletStart + R - cos22half * IR, sin22half * IR, highZ],
        )
        inletTopBlock, inletTopPoints = mesh.custom(
            [
                [inletStart + R - iScale * IR, iScale * IR, lowZ],
                [inletStart + R + iScale * IR, iScale * IR, lowZ],
                [inletStart + R + cos45 * IR, cos45 * IR, lowZ],
                [inletStart + R - cos45 * IR, cos45 * IR, lowZ],
                [inletStart + R - iScale * IR, iScale * IR, highZ],
                [inletStart + R + iScale * IR, iScale * IR, highZ],
                [inletStart + R + cos45 * IR, cos45 * IR, highZ],
                [inletStart + R - cos45 * IR, cos45 * IR, highZ],
            ],
            ["membrane", "symm", "chamberWall", "chamberWall", "chamberWall", "inlet"],
            [fine * 12, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(inletTopPoints[2], inletTopPoints[3], [inletStart + R, IR, lowZ])
        mesh.arc(inletTopPoints[6], inletTopPoints[7], [inletStart + R, IR, highZ])
        inletRightBlock, inletRightPoints = mesh.custom(
            [
                [inletStart + R + iScale * IR, 0, lowZ],
                [inletStart + R + IR, 0, lowZ],
                [inletStart + R + cos45 * IR, cos45 * IR, lowZ],
                [inletStart + R + iScale * IR, iScale * IR, lowZ],
                [inletStart + R + iScale * IR, 0, highZ],
                [inletStart + R + IR, 0, highZ],
                [inletStart + R + cos45 * IR, cos45 * IR, highZ],
                [inletStart + R + iScale * IR, iScale * IR, highZ],
            ],
            ["membrane", "symm", "chamberWall", "chamberWall", "chamberWall", "inlet"],
            [fine * 5, fine * 5, fine * cells],
            [1, symmScale, scaling],
        )
        mesh.arc(
            inletRightPoints[1],
            inletRightPoints[2],
            [inletStart + R + cos22half * IR, sin22half * IR, lowZ],
        )
        mesh.arc(
            inletRightPoints[5],
            inletRightPoints[6],
            [inletStart + R + cos22half * IR, sin22half * IR, highZ],
        )


def generate(
    path,
    *,
    flow=50,
    length=0.085,
    width=0.039,
    height=0.00225,
    diameter=0.0055,
    corner_radius=0.003,
    offset=0.0045,
    distributor_width=0.0295,
    chamber_height=0.01925,
    refinement=1,
    membrane_expansion=5,
    draw=0.065,
    feed=0.00065,
    resistance=150666,
    end=60000,
):
    """Write the source-derived quarter mesh and canonical two-sided FO fields."""
    import json
    from pathlib import Path
    from .channel import generate as channel
    from .chambers import configure_2016

    values = (
        flow,
        length,
        width,
        height,
        diameter,
        corner_radius,
        offset,
        distributor_width,
        chamber_height,
        refinement,
        membrane_expansion,
        draw,
        feed,
        resistance,
        end,
    )
    if not all(math.isfinite(v) for v in values):
        raise ValueError("All geometry, mesh and operating parameters must be finite")
    inlet_height = chamber_height - height
    if not (
        flow > 0
        and 0 <= feed < draw <= 0.09
        and resistance >= 0
        and end >= 1
        and refinement >= 0.25
        and membrane_expansion >= 1
        and height > 0
        and inlet_height > height
        and 0 < diameter / 2 < corner_radius
        and 0 < offset < length / 2 - 2 * corner_radius - 0.0005
        and corner_radius + diameter / 2 < distributor_width / 2
        and width / 2 - corner_radius - distributor_width / 2 >= 1 / 3000
    ):
        raise ValueError(
            "Dimensions must leave room for the original CF042 blocks and ports"
        )
    mesh = Geometry(
        length=length,
        width=width,
        height=height,
        diameter=diameter,
        corner_radius=corner_radius,
        offset=offset,
        distributor_width=distributor_width,
        inlet_height=inlet_height,
        refinement=refinement,
        membrane_expansion=membrane_expansion,
    ).mesh
    path = Path(path)
    channel(path, end=end)
    (path / "system/blockMeshDict").write_text(mesh.dictionary())
    configure_2016(path, flow, draw, feed, resistance, end)
    layers = max(1, round(25 * refinement))
    log_step = math.log(membrane_expansion) / (layers - 1)
    first_layer = (
        height / layers
        if membrane_expansion == 1
        else height
        * math.expm1(log_step)
        * math.exp(-layers * log_step)
        / -math.expm1(-layers * log_step)
    )
    (path / "case.json").write_text(
        json.dumps(
            dict(
                geometry="CF042 toolkit geometry",
                paper="10.1016/j.seppur.2015.12.017",
                geometry_source="MembraneSimKit/blockCF042",
                geometry_source_revision="2c4344144dadc3f2b486432724468a91652bbe3a",
                flow_ml_min=flow,
                length_m=length,
                width_m=width,
                height_m=height,
                inlet_offset_m=offset,
                inlet_diameter_m=diameter,
                chamber_height_m=chamber_height,
                channel_corner_radius_m=corner_radius,
                distributor_width_m=distributor_width,
                inlet_height_m=inlet_height,
                inlet_buffer_distance_m=0.0005,
                mesh_merge_gap_m=0.000001,
                draw_mass_fraction=draw,
                feed_mass_fraction=feed,
                support_resistance_s_m=resistance,
                refinement=refinement,
                half_width=True,
                wall_normal_expansion=membrane_expansion,
                membrane_layers=layers,
                nominal_membrane_cell_height_m=first_layer,
            ),
            indent=2,
        )
        + "\n"
    )


def prepare_mesh(path, length):
    """Merge the quarter first, then mirror and pair its membrane faces."""
    import json
    import shutil
    import subprocess
    from pathlib import Path

    path = Path(path).resolve()
    metadata = json.loads((path / "case.json").read_text())
    extent = max(length, metadata["width_m"], metadata["chamber_height_m"]) + 1
    initial = path / "initialFields"
    (path / "0").rename(initial)

    def tool(name, *options):
        with (path / ("log." + name)).open("a") as log:
            subprocess.run(
                [name, "-case", str(path), *options],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )

    tool("blockMesh")
    for point, normal in [(f"{length / 2} 0 0", "1 0 0"), ("0 0 0", "0 0 1")]:
        (path / "system/mirrorMeshDict").write_text(header("mirrorMeshDict") + f"""
planeType pointAndNormal;
pointAndNormalDict {{ basePoint ({point}); normalVector ({normal}); }}
planeTolerance 1e-8;
""")
        tool("mirrorMesh", "-overwrite")
    ports = [
        ("feedIn", -extent, length / 2, -extent, 0),
        ("feedOut", length / 2, extent, -extent, 0),
        ("drawIn", length / 2, extent, 0, extent),
        ("drawOut", -extent, length / 2, 0, extent),
    ]
    actions = [
        f"{{ name membraneFaces; type faceSet; action new; source boxToFace; box ({-extent} {-extent} -1e-10) ({extent} {extent} 1e-10); }}",
        "{ name membraneZone; type faceZoneSet; action new; source setToFaceZone; faceSet membraneFaces; }",
    ]
    for name, lo, hi, zlo, zhi in ports:
        actions += [
            f"{{ name {name}; type faceSet; action new; source patchToFace; patches (inlet); }}",
            f"{{ name {name}; type faceSet; action subset; source boxToFace; box ({lo} {-extent} {zlo}) ({hi} {extent} {zhi}); }}",
        ]
    (path / "system/topoSetDict").write_text(
        header("topoSetDict") + "actions (" + "\n".join(actions) + ");\n"
    )
    tool("topoSet")
    (path / "system/createBafflesDict").write_text(header("createBafflesDict") + """
internalFacesOnly true;
baffles { membrane { type faceZone; zoneName membraneZone;
 patches { master { name membrane; type wall; } slave { name membrane; type wall; } }
} }
""")
    tool("createBaffles", "-overwrite")
    tool("topoSet")
    patches = [
        f"{{ name {name}; patchInfo {{ type patch; }} constructFrom set; set {name}; }}"
        for name, *_ in ports
    ]
    patches += [
        "{ name walls; patchInfo { type wall; } constructFrom patches; patches (chamberWall inletWall tempWall1 tempWall2); }",
        "{ name sides; patchInfo { type symmetryPlane; } constructFrom patches; patches (symm); }",
    ]
    (path / "system/createPatchDict").write_text(
        header("createPatchDict")
        + "pointSync false; patches ("
        + "\n".join(patches)
        + ");\n"
    )
    tool("createPatch", "-overwrite")
    tool("checkMesh")
    if "Mesh OK" not in (path / "log.checkMesh").read_text():
        raise RuntimeError("Original CF042 mesh check failed")
    (path / "0").mkdir(exist_ok=True)
    for field in initial.iterdir():
        shutil.copy2(field, path / "0" / field.name)
    shutil.rmtree(initial)
