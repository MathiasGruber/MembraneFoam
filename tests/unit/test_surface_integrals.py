"""Area weighting on a deliberately nonuniform membrane mesh."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.results import surface_integrals


class SurfaceIntegrationTests(unittest.TestCase):
    def test_nonuniform_faces(self):
        with tempfile.TemporaryDirectory() as directory:
            surface = Path(directory) / "surface.csv"
            surface.write_text(
                "nz,mass_fraction,area_m2\n-1,0.02,3\n-1,0.06,1\n1,0.001,4\n"
            )
            result = surface_integrals(surface)
        self.assertAlmostEqual(result["negative_z"]["area_m2"], 4)
        self.assertAlmostEqual(result["negative_z"]["mass_fraction_integral_m2"], 0.12)
        self.assertAlmostEqual(result["negative_z"]["mean_mass_fraction"], 0.03)
        self.assertAlmostEqual(result["positive_z"]["mean_mass_fraction"], 0.001)
