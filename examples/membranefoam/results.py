#!/usr/bin/env python3
"""Produce a machine-readable run record from solver logs; reject incomplete runs."""

import csv
import hashlib
import json
import math
import re
from pathlib import Path


def summarize(case):
    case = Path(case)
    log = (case / "log.solver").read_text()
    if not re.search(r"^End\s*$", log, re.M) or "FOAM FATAL" in log:
        raise ValueError("Solver did not complete successfully")
    # Older solvers can divide by zero in this wall-clock performance report.
    numerical_log = re.sub(r"^Simulation Speed = .*?$", "", log, flags=re.M)
    if re.search(r"(?<![A-Za-z])[-+]?(?:nan|inf)(?![A-Za-z])", numerical_log, re.I):
        raise ValueError("Non-finite value in solver log")

    def last(pattern):
        matches = re.findall(pattern, log, re.M)
        return float(matches[-1]) if matches else None

    fluxes = [float(v) for v in re.findall(r"Water flux, \w+: ([\d.eE+-]+)", log)]
    window = fluxes[-100:]
    record = json.loads((case / "case.json").read_text())
    record.update(
        final_time=last(r"^Time = ([\d.eE+-]+)"),
        solute_max=last(r"m_A max/min : ([\d.eE+-]+)"),
        solute_min=last(r"m_A max/min : [\d.eE+-]+ ([\d.eE+-]+)"),
        flux_window_samples=len(window),
        water_mass_flux_kg_m2_h=last(r"Water flux, \w+: ([\d.eE+-]+)"),
        membrane_salt_imbalance_kg_s=last(
            r"patch=membrane mass_kg_s=[\d.eE+-]+ salt_kg_s=([\d.eE+-]+)"
        ),
        salt_flux_g_m2_h=last(r"Salt Flux: ([\d.eE+-]+)"),
        mass_imbalance_kg_s=last(r"total_mass_kg_s=([\d.eE+-]+)"),
        salt_imbalance_kg_s=last(
            r"MEMBRANE_CHECK total_mass_kg_s=[\d.eE+-]+ total_salt_kg_s=([\d.eE+-]+)"
        ),
        discrete_salt_imbalance_kg_s=last(
            r"MEMBRANE_DISCRETE total_salt_kg_s=([\d.eE+-]+)"
        ),
        membrane_salt_transfer_kg_s=last(r"membrane_salt_transfer_kg_s=([\d.eE+-]+)"),
        relative_salt_imbalance=last(r"relative_salt_imbalance=([\d.eE+-]+)"),
        relative_mass_imbalance=last(r"relative_mass_imbalance=([\d.eE+-]+)"),
        converged="SIMPLE solution converged" in log,
        last_100_flux_relative_span=(
            (max(window) - min(window)) / max(abs(window[-1]), 1e-30)
            if window
            else None
        ),
        openfoam_build=re.search(r"^Build\s*:\s*(.+)", log, re.M).group(1),
    )
    if record["solute_min"] is None:
        record["solute_min"] = last(r"^m_A: min = ([\d.eE+-]+)")
        record["solute_max"] = last(r"^m_A: min = [\d.eE+-]+\s+max = ([\d.eE+-]+)")
    if "MEMBRANE_FLUX" in log:
        record["water_mass_flux_kg_m2_h"] = last(r"water_kg_m2_h=([\d.eE+-]+)")
        record["salt_flux_g_m2_h"] = last(r"salt_g_m2_h=([\d.eE+-]+)")
    times = case / "log.time"
    if times.exists():
        m = re.search(r"^real ([\d.]+)", times.read_text(), re.M)
        if m:
            record["wall_seconds"] = float(m.group(1))
    for field in ("p", "m_A", "Ux", "Uy", "Uz"):
        record[field + "_initial_residual"] = last(
            r"Solving for " + field + r", Initial residual = ([\d.eE+-]+)"
        )
    history = transient_history(log)
    record["run_type"] = "transient" if history else "steady"
    if record.get("end_time_s") is not None and not history:
        raise ValueError("Transient run is missing accumulation-inclusive audits")
    if history:
        # Assess the physical trajectory, using the final PISO correction at
        # each time step rather than intermediate nonlinear iterates.
        bounds, current = [], None
        for match in re.finditer(
            r"^Time = .*?$|^m_A: min = ([\d.eE+-]+)\s+max = ([\d.eE+-]+)",
            log, re.M,
        ):
            if match.group(1) is None:
                if current is not None:
                    bounds.append(current)
                current = None
            else:
                current = tuple(float(v) for v in match.groups())
        if current is not None:
            bounds.append(current)
        if not bounds:
            raise ValueError("Transient run is missing time-step concentration bounds")
        record["transient_solute_min"] = min(row[0] for row in bounds)
        record["transient_solute_max"] = max(row[1] for row in bounds)
        record["transient_audit"] = dict(
            samples=len(history),
            first_time_s=history[0]["time"],
            last_time_s=history[-1]["time"],
            max_relative_salt_balance=max(
                row["relative_salt_balance"] for row in history
            ),
            max_relative_mass_balance=max(
                row["relative_mass_balance"] for row in history
            ),
        )
        record["converged"] = False
    return record


def transient_history(log):
    """Extract write-time balances and the final FO flux preceding each audit."""
    history = []
    flux = None
    for line in log.splitlines():
        match = re.search(r"Water flux, \w+: ([\d.eE+-]+)", line)
        if match:
            flux = float(match[1])
        if line.startswith("MEMBRANE_TRANSIENT "):
            row = {
                key: float(value) for key, value in re.findall(r"(\w+)=([^\s]+)", line)
            }
            if flux is not None:
                row["water_kg_m2_h"] = flux
            if not all(math.isfinite(value) for value in row.values()):
                raise ValueError("Non-finite transient audit")
            history.append(row)
    return history


def surface_integrals(surface):
    """Integrate membrane concentration, grouping faces by their z-normal sign."""
    groups = {}
    with Path(surface).open() as source:
        rows = csv.DictReader(source)
        if "area_m2" not in (rows.fieldnames or []):
            return None  # Older exports do not contain integration weights.
        for row in rows:
            area, value, nz = (
                float(row[k]) for k in ("area_m2", "mass_fraction", "nz")
            )
            if not all(math.isfinite(v) for v in (area, value, nz)) or area <= 0:
                raise ValueError("Invalid membrane surface area or concentration")
            side = (
                "positive_z"
                if nz > 1e-12
                else "negative_z" if nz < -1e-12 else "parallel_z"
            )
            group = groups.setdefault(
                side, dict(area_m2=0.0, mass_fraction_integral_m2=0.0)
            )
            group["area_m2"] += area
            group["mass_fraction_integral_m2"] += area * value
    for group in groups.values():
        group["mean_mass_fraction"] = (
            group["mass_fraction_integral_m2"] / group["area_m2"]
        )
    return groups


def collect(case, output):
    """Export a converged steady run or a completed, audited transient run."""
    case, output = Path(case), Path(output)
    record = summarize(case)
    if (
        record["solute_min"] is None
        or record["solute_max"] is None
        or record["solute_min"] < -1e-6
        or record["solute_max"] > 1
    ):
        raise ValueError("Unphysical salt mass fraction; do not publish this run")
    if record["run_type"] == "transient":
        if record["transient_solute_min"] < -1e-6 or record["transient_solute_max"] > 1:
            raise ValueError("Unphysical salt mass fraction in transient trajectory")
        audit = record["transient_audit"]
        if (
            audit["max_relative_salt_balance"] >= 0.001
            or audit["max_relative_mass_balance"] >= 1e-6
        ):
            raise ValueError(
                "Transient accumulation-inclusive conservation check failed"
            )
        if (
            audit["last_time_s"]
            < record.get("end_time_s", audit["last_time_s"]) - 1e-10
        ):
            raise ValueError("Transient run did not reach its requested end time")
    else:
        if not record["converged"]:
            raise ValueError(
                "SIMPLE did not converge; do not publish this run as a reproduced result"
            )
        if (
            record["relative_salt_imbalance"] is None
            or record["relative_salt_imbalance"] >= 0.01
        ):
            raise ValueError(
                "Salt imbalance must be below 1% of membrane transfer; tighten convergence"
            )
        if (
            record["relative_mass_imbalance"] is None
            or record["relative_mass_imbalance"] >= 1e-6
        ):
            raise ValueError("Relative mass imbalance must be below 1e-6")
    record["configuration"] = {
        str(path.relative_to(case)): path.read_text()
        for path in (
            case / "system/controlDict",
            case / "system/fvSolution",
            case / "system/fvSchemes",
            case / "constant/transportProperties",
        )
    }
    mesh_log = case / "log.checkMesh"
    if mesh_log.exists():
        mesh_text = mesh_log.read_text()
        record["mesh_check_passed"] = "Mesh OK" in mesh_text
        count = re.search(r"^\s*cells:\s+(\d+)", mesh_text, re.M)
        if count:
            record["cells"] = int(count.group(1))
    record["solver_log_sha256"] = hashlib.sha256(
        (case / "log.solver").read_bytes()
    ).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    if record["run_type"] == "transient":
        history = transient_history((case / "log.solver").read_text())
        with output.with_name(output.name + "-history.csv").open("w") as target:
            writer = csv.DictWriter(
                target, fieldnames=list(history[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(history)
    parts = sorted(case.glob("processor[0-9]*/membrane.csv"))
    if not parts:
        parts = [case / "membrane.csv"]
    with output.with_name(output.name + "-surface.csv").open("w") as target:
        writer = csv.writer(target, lineterminator="\n")
        header = None
        for part in parts:
            with part.open() as source:
                rows = csv.reader(source)
                columns = next(rows)
                if header is None:
                    header = columns
                    writer.writerow(columns)
                elif columns != header:
                    raise ValueError(f"Inconsistent columns: {part}")
                writer.writerows(rows)
    integrals = surface_integrals(output.with_name(output.name + "-surface.csv"))
    if integrals is not None:
        record["surface_integrals"] = integrals
    output.with_name(output.name + ".json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n"
    )
    return record


def refinement(root):
    """Estimate observed order and GCI from channel-40/80/160 records."""
    root = Path(root)
    records = [
        json.loads((root / f"channel-{n}.json").read_text()) for n in (40, 80, 160)
    ]
    if not all(r["converged"] for r in records):
        raise ValueError("All three runs must converge")
    c, m, f = [r["water_mass_flux_kg_m2_h"] for r in records]
    if (c - m) * (m - f) <= 0:
        raise ValueError("Non-monotone sequence: do not apply this GCI estimate")
    p = math.log(abs((c - m) / (m - f))) / math.log(2)
    if p <= 0:
        raise ValueError("No positive observed order")
    result = dict(
        geometry="2D verification channel only",
        refinement_ratio=2,
        observed_order=p,
        fine_grid_gci_percent=100 * 1.25 * abs((m - f) / f) / (2**p - 1),
        coarse_to_medium_change_percent=100 * abs((c - m) / m),
        medium_to_fine_change_percent=100 * abs((m - f) / f),
        extrapolated_water_flux_kg_m2_h=f + (f - m) / (2**p - 1),
        safety_factor=1.25,
    )
    (root / "channel-refinement.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
