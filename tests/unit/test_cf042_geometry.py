"""The maintained construction must preserve the qualified source quarter mesh."""

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.cf042_original import Geometry, generate


class OriginalCF042Tests(unittest.TestCase):
    def test_source_geometry_fingerprint(self):
        mesh = Geometry(
            length=0.085,
            width=0.039,
            height=0.00225,
            diameter=0.0055,
            corner_radius=0.003,
            offset=0.0045,
            distributor_width=0.0295,
            inlet_height=0.017,
            refinement=1,
        ).mesh
        # Fingerprint of the independently translated MembraneSimKit quarter
        # mesh that passed v2606 checkMesh before integration into this runner.
        self.assertEqual(
            hashlib.sha256(mesh.dictionary().encode()).hexdigest(),
            "124ee1cc22117d35b1bfe556ed050da4e6904796e42b96026b40fd797f4585d4",
        )
        self.assertEqual(sum(n[0] * n[1] * n[2] for _, n, _ in mesh.blocks), 206570)

    def test_membrane_grading_changes_resolution_without_moving_geometry(self):
        parameters = dict(
            length=0.085,
            width=0.039,
            height=0.00225,
            diameter=0.0055,
            corner_radius=0.003,
            offset=0.0045,
            distributor_width=0.0295,
            inlet_height=0.017,
            refinement=1,
        )
        source = Geometry(**parameters).mesh
        graded = Geometry(**parameters, membrane_expansion=100).mesh
        self.assertEqual(source.points, graded.points)
        self.assertEqual([b[:2] for b in source.blocks], [b[:2] for b in graded.blocks])
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory) / "case"
            generate(case, membrane_expansion=100)
            record = json.loads((case / "case.json").read_text())
            self.assertEqual(record["membrane_layers"], 25)
            self.assertGreater(record["nominal_membrane_cell_height_m"], 3e-6)
            self.assertLess(record["nominal_membrane_cell_height_m"], 5e-6)

    def test_original_inputs_are_recorded_and_invalid_ports_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory) / "case"
            generate(case)
            record = json.loads((case / "case.json").read_text())
            self.assertEqual(record["geometry"], "CF042 toolkit geometry")
            self.assertEqual(record["channel_corner_radius_m"], 0.003)
            self.assertEqual(record["distributor_width_m"], 0.0295)
            self.assertAlmostEqual(record["inlet_height_m"], 0.017)
            with self.assertRaises(ValueError):
                generate(Path(directory) / "bad", diameter=0.008)
            self.assertFalse((Path(directory) / "bad").exists())
