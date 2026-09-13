"""The public coupling controls must configure both the solver and its record."""

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

RUNNER = Path(__file__).resolve().parents[2] / "examples/run.py"


class TransientGenerationTests(unittest.TestCase):
    def test_default_and_source_coupling(self):
        for options, correctors, predictor in [
            ([], 2, True),
            (["--pressure-correctors", "1", "--no-momentum-predictor"], 1, False),
        ]:
            with self.subTest(options=options), tempfile.TemporaryDirectory() as folder:
                case = Path(folder) / "case"
                subprocess.run([
                    sys.executable, str(RUNNER), "channel", str(case),
                    "--generate-only", "--transient", "0.001", "--time-step", "0.0001",
                    *options,
                ], check=True, capture_output=True, text=True)
                config = (case / "system/fvSolution").read_text()
                piso = re.search(r"PISO\s*\{([^}]+)\}", config)[1]
                self.assertRegex(piso, rf"\bnCorrectors\s+{correctors};")
                self.assertRegex(piso, rf"\bmomentumPredictor\s+{str(predictor).lower()};")
                record = json.loads((case / "case.json").read_text())
                self.assertEqual(record["n_pressure_correctors"], correctors)
                self.assertIs(record["momentum_predictor"], predictor)

    def test_invalid_or_steady_only_usage(self):
        for options in [
            ["--no-momentum-predictor"],
            ["--pressure-correctors", "1"],
            ["--transient", "0.001", "--time-step", "0.0001", "--pressure-correctors", "0"],
        ]:
            with self.subTest(options=options), tempfile.TemporaryDirectory() as folder:
                case = Path(folder) / "case"
                result = subprocess.run([
                    sys.executable, str(RUNNER), "channel", str(case), "--generate-only",
                    *options,
                ], capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertFalse(case.exists())


if __name__ == "__main__":
    unittest.main()
