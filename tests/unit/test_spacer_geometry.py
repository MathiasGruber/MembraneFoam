"""Check spacer displacement from generated straight block geometry."""

import json
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.block import generate
from membranefoam.spacers import wall_blocks


class SpacerGeometryTests(unittest.TestCase):
    def test_triangle_fluid_volume_matches_requested_displacement(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory) / "triangle"
            generate(case, spacers=2, spacer_shape="triangle", spacer_area_fraction=0.2)
            text = (case / "system/blockMeshDict").read_text()
            coordinates = text.split("vertices (\n", 1)[1].split("\n);", 1)[0]
            points = [
                tuple(map(float, row.strip("()").split()))
                for row in coordinates.splitlines()
            ]
            volume = 0.0
            for row in re.findall(r"hex \(([^)]+)\)", text):
                vertices = [points[int(index)] for index in row.split()]
                if min(p[0] for p in vertices) < 0.002 - 1e-12:
                    continue
                # All main-chamber blocks are straight extrusions in y.
                polygon = [vertices[i] for i in (0, 1, 5, 4)]
                area = (
                    abs(
                        sum(
                            a[0] * b[2] - b[0] * a[2]
                            for a, b in zip(polygon, polygon[1:] + polygon[:1])
                        )
                    )
                    / 2
                )
                volume += area * (
                    max(p[1] for p in vertices) - min(p[1] for p in vertices)
                )
            expected = (0.04 - 0.002) * 0.02 * 0.002 - 0.2 * 0.002**2 * 0.02
            self.assertAlmostEqual(volume, expected, delta=1e-16)
            record = json.loads((case / "case.json").read_text())
            self.assertEqual(record["spacer_shape"], "triangle")
            self.assertAlmostEqual(
                record["spacer_triangle_side_m"] ** 2 * math.sqrt(3) / 4, 0.2 * 0.002**2
            )

    def test_equal_area_cylinder_uses_the_same_cross_section(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory) / "cylinder"
            generate(case, spacers=2, spacer_area_fraction=0.2)
            record = json.loads((case / "case.json").read_text())
            self.assertAlmostEqual(
                math.pi * record["spacer_radius_m"] ** 2, 0.2 * 0.002**2
            )

    def test_wall_boundaries_displace_the_target_area(self):
        for shape in ("cylinder", "triangle"):
            for placement in ("up", "down"):
                for fraction in (0.1, 0.2, 0.3):
                    with self.subTest(
                        shape=shape, placement=placement, fraction=fraction
                    ):
                        profiles, footprint = wall_blocks(
                            shape,
                            placement,
                            fraction,
                            [0, 0.00045, 0.001, 0.002],
                            lambda lo, hi: (7, 2),
                        )
                        radius = math.sqrt(fraction * 0.002**2 / (0.75 * math.pi + 0.5))
                        center_z = (
                            radius / math.sqrt(2)
                            if placement == "down"
                            else 0.002 - radius / math.sqrt(2)
                        )
                        area = membrane_width = 0
                        for profile in profiles:
                            q = profile["quad"]
                            curved = {
                                tuple(sorted((a, b))) for a, b, _ in profile["arcs"]
                            }
                            cell_area = 0
                            for i, j in zip(range(4), (1, 2, 3, 0)):
                                a, b = q[i], q[j]
                                if tuple(sorted((i, j))) in curved:
                                    start = math.atan2(a[1] - center_z, a[0])
                                    end = math.atan2(b[1] - center_z, b[0])
                                    delta = (end - start + math.pi) % (
                                        2 * math.pi
                                    ) - math.pi
                                    cell_area += (
                                        -center_z * (b[0] - a[0]) + radius**2 * delta
                                    ) / 2
                                else:
                                    cell_area += (a[0] * b[1] - b[0] * a[1]) / 2
                                if abs(a[1]) < 1e-12 and abs(b[1]) < 1e-12:
                                    membrane_width += abs(b[0] - a[0])
                            self.assertGreater(cell_area, 0)
                            area += cell_area
                        self.assertAlmostEqual(
                            area, (1 - fraction) * 0.002**2, delta=1e-18
                        )
                        expected_width = 0.002 - (
                            footprint if placement == "down" else 0
                        )
                        self.assertAlmostEqual(
                            membrane_width, expected_width, delta=1e-15
                        )

    def test_wall_spacer_area_and_source_changing_order(self):
        with tempfile.TemporaryDirectory() as directory:
            for shape in ("cylinder", "triangle"):
                case = Path(directory) / shape
                generate(
                    case,
                    spacers=6,
                    spacer_shape=shape,
                    spacer_area_fraction=0.2,
                    spacer_placement="changing",
                )
                record = json.loads((case / "case.json").read_text())
                self.assertEqual(
                    record["spacer_placements"],
                    ["down", "up", "down", "down", "up", "down"],
                )
                self.assertEqual(record["spacer_membrane_contact_count"], 4)
                self.assertAlmostEqual(record["spacer_cross_section_m2"], 0.8e-6)
                footprint = (
                    math.sqrt(2 * 0.8e-6 / (0.75 * math.pi + 0.5))
                    if shape == "cylinder"
                    else math.sqrt(4 * 0.8e-6 / math.sqrt(3))
                )
                self.assertAlmostEqual(
                    record["expected_active_membrane_area_full_m2"],
                    (0.08 - 4 * footprint) * 0.04,
                )

    def test_invalid_shape_inputs_do_not_create_a_case(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory) / "invalid"
            for parameters in [
                dict(spacers=2, spacer_shape="triangle"),
                dict(
                    spacers=2,
                    spacer_shape="triangle",
                    spacer_radius=0.0005,
                    spacer_area_fraction=0.2,
                ),
                dict(spacers=2, spacer_area_fraction=float("nan")),
                dict(spacer_area_fraction=0.2),
            ]:
                with self.assertRaises(ValueError):
                    generate(case, **parameters)
                self.assertFalse(case.exists())


if __name__ == "__main__":
    unittest.main()
