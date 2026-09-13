"""Known convergence laws and unsuitable real sequences for conditional GCI."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.results import effective_grid_estimate


class GridEstimateTests(unittest.TestCase):
    def test_nasa_pressure_recovery_example(self):
        result = effective_grid_estimate([1, 8, 64], [0.96178, 0.96854, 0.97050])
        self.assertAlmostEqual(result["observed_order"], 1.786170, places=6)
        self.assertAlmostEqual(result["extrapolated_value"], 0.971300, places=6)

    def test_known_nonuniform_refinement_law(self):
        counts = [1000, 3375, 8000]
        values = [7 + 2 * (counts[0] / n) ** (2.3 / 3) for n in counts]
        result = effective_grid_estimate(counts, values)
        self.assertAlmostEqual(result["observed_order"], 2.3, places=11)
        self.assertAlmostEqual(result["extrapolated_value"], 7, places=12)
        expected = 100 * 1.25 * abs(values[-1] - 7) / values[-1]
        self.assertAlmostEqual(result["fine_grid_gci_percent"], expected, places=11)

    def test_unsuitable_observed_sequences(self):
        sequences = [
            ([1066000, 2288976, 3698800],
             [11.9738916926, 11.9734582896, 11.9681438398]),
            ([1140744, 1926360, 4428480], [5.62717646, 5.62805830, 5.62574605]),
            ([1000, 3375, 8000], [1, 1, 1]),
        ]
        for counts, values in sequences:
            with self.subTest(values=values), self.assertRaises(ValueError):
                effective_grid_estimate(counts, values)

    def test_invalid_inputs(self):
        for counts, values in [
            ([1, 1, 64], [3, 2, 1]),
            ([1, 8], [3, 2, 1]),
            ([1, 8, 64], [3, float("nan"), 1]),
            ([1, 8, 64], [3, 2, 0]),
            ([1, 8, 64], [3, 0, -1]),
            ([1, 8.5, 64], [3, 2, 1]),
            ([100, 101, 102], [3, 2, 1]),
        ]:
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                effective_grid_estimate(counts, values)


if __name__ == "__main__":
    unittest.main()
