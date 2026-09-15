"""Independent checks against the 2016 analytical figure, including its units."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.paper_models import film_coefficient, water_flux_mass_fractions


class PaperModels(unittest.TestCase):
    def test_equal_concentrations(self):
        self.assertEqual(water_flux_mass_fractions(0, 1.6e-12, 0, 150666), 0)
        self.assertEqual(
            water_flux_mass_fractions(0.065, 1.6e-12, 8e-8, 150666, 0.065), 0
        )

    def test_archived_analytical_curve(self):
        # Figure 2 analytical curves, read from the author figure PDF. The 0.06
        # tolerance covers raster resolution and dashed-line sampling in kg/m²/h.
        for flow, length, expected in [
            (50, 0.085, 10.9722),
            (50, 0.2, 10.0833),
            (50, 0.5, 9.0556),
            (5, 0.085, 8.4306),
            (5, 0.2, 7.4444),
            (500, 0.2, 12.2778),
        ]:
            with self.subTest(flow=flow, length=length):
                j = water_flux_mass_fractions(
                    0.065,
                    1.61111e-12,
                    8.33333e-8,
                    150666,
                    0.00065,
                    film_coefficient(flow, length, 0.039, 0.00225),
                )
                # The publication's analytical plot uses reference density 1000.
                self.assertAlmostEqual(j * 1000 * 3600, expected, delta=0.06)

    def test_invalid_concentrations(self):
        for feed, draw in [(-0.01, 0.065), (0.07, 0.065), (0, 0.1)]:
            with self.assertRaises(ValueError):
                water_flux_mass_fractions(draw, 1.6e-12, 8e-8, 150666, feed)


if __name__ == "__main__":
    unittest.main()
