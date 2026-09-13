"""Transient samples must use the final update, not duplicate PISO correctors."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
from membranefoam.results import collect, summarize, transient_history


class TransientHistoryTests(unittest.TestCase):
    def test_final_flux_at_each_audit(self):
        log = """
Time = 0.001
    Water flux, advanced: 15 kg/(h*m2)
    Water flux, advanced: 12 kg/(h*m2)
MEMBRANE_TRANSIENT time=0.001 relative_salt_balance=1e-5 relative_mass_balance=1e-12
Time = 0.002
    Water flux, advanced: 11 kg/(h*m2)
MEMBRANE_TRANSIENT time=0.002 relative_salt_balance=2e-5 relative_mass_balance=2e-12
"""
        rows = transient_history(log)
        self.assertEqual([row["water_kg_m2_h"] for row in rows], [12, 11])
        self.assertEqual([row["time"] for row in rows], [0.001, 0.002])

    def test_nonfinite_audit_rejected(self):
        with self.assertRaises(ValueError):
            transient_history("MEMBRANE_TRANSIENT time=0.001 relative_salt_balance=nan")

    def test_piso_summary_and_legacy_timing(self):
        log = """Build : test
Time = 0.002
m_A: min = 0    max = 0.02873
    Water flux, advanced: 12 kg/(h*m2)
MEMBRANE_TRANSIENT time=0.002 relative_salt_balance=2e-5 relative_mass_balance=2e-12
Simulation Speed = inf s / day
End
"""
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            (case / "case.json").write_text('{"end_time_s": 0.002}')
            (case / "log.solver").write_text(log)
            result = summarize(case)
            self.assertEqual(result["run_type"], "transient")
            self.assertFalse(result["converged"])
            self.assertEqual(result["solute_min"], 0)
            self.assertEqual(result["solute_max"], 0.02873)
            (case / "log.solver").write_text(log.replace("min = 0", "min = nan"))
            with self.assertRaises(ValueError):
                summarize(case)

    def test_decimal_result_stems_remain_distinct(self):
        log = """Build : test
Time = 0.002
m_A: min = 0    max = 0.02873
    Water flux, advanced: 12 kg/(h*m2)
MEMBRANE_TRANSIENT time=0.002 relative_salt_balance=2e-5 relative_mass_balance=2e-12
End
"""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = root / "case"
            case.mkdir()
            (case / "case.json").write_text('{"end_time_s": 0.002}')
            (case / "log.solver").write_text(log)
            (case / "membrane.csv").write_text(
                "nz,mass_fraction,area_m2\n-1,0.02,1\n1,0.001,1\n"
            )
            for name in (
                "system/controlDict", "system/fvSolution", "system/fvSchemes",
                "constant/transportProperties",
            ):
                path = case / name
                path.parent.mkdir(exist_ok=True)
                path.write_text("// test dictionary\n")
            for stem in ("step-0.0005", "step-0.001"):
                collect(case, root / stem)
            records = sorted(root.glob("*.json"))
            self.assertEqual(
                [p.name for p in records], ["step-0.0005.json", "step-0.001.json"]
            )
            for record in records:
                for suffix in ("-history.csv", "-surface.csv"):
                    self.assertTrue(record.with_name(record.stem + suffix).is_file())

    def test_unphysical_earlier_step_cannot_be_hidden_by_recovery(self):
        log = """Build : test
Time = 0.001
m_A: min = -0.01    max = 0.07
Time = 0.002
m_A: min = 0    max = 0.07
MEMBRANE_TRANSIENT time=0.002 relative_salt_balance=1e-5 relative_mass_balance=1e-12
End
"""
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            (case / "case.json").write_text('{"end_time_s": 0.002}')
            (case / "log.solver").write_text(log)
            self.assertEqual(summarize(case)["solute_min"], 0)
            with self.assertRaisesRegex(ValueError, "transient trajectory"):
                collect(case, case / "result")
            # An intermediate correction is not a completed physical time step.
            (case / "log.solver").write_text(log.replace(
                "Time = 0.002", "m_A: min = 0    max = 0.07\nTime = 0.002"
            ))
            self.assertEqual(summarize(case)["transient_solute_min"], 0)
