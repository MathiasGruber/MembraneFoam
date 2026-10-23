#!/usr/bin/env python3
"""Render reproducible publication comparisons. Requires numpy and matplotlib."""

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .paper_models import film_coefficient, water_flux, mass_fraction
from .results import paired_face_indices

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/figures"
RESULTS = ROOT / "runs/results"
INK = "#16324a"
TEAL = "#087f8c"
ORANGE = "#dc7542"
MUTED = "#677c89"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "svg.hashsalt": "MembraneFoam-v2606",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": INK,
        "text.color": INK,
        "axes.edgecolor": "#bdcbd3",
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.titleweight": "bold",
        "figure.facecolor": "#fafcfb",
        "axes.facecolor": "#fafcfb",
        "savefig.facecolor": "#fafcfb",
        "grid.color": "#dce5e9",
        "grid.alpha": 0.65,
    }
)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / (name + ".svg"), bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def analytical():
    """Compare archived CFD, the published film model, and validated CF042 runs."""
    from matplotlib.lines import Line2D
    from .paper_models import PAPER_2016_DIMENSIONS, water_flux_mass_fractions

    # Only converged, conserved runs with the reference operating point are eligible.
    candidates = {}
    for source in RESULTS.glob("cf042-*.json"):
        record = json.loads(source.read_text())
        if (
            not record.get("converged")
            or record.get("geometry", "CF042 reconstruction")
            not in ("CF042 reconstruction", "CF042 toolkit geometry")
            or record.get("relative_salt_imbalance", 1) >= 0.001
            or record.get("relative_mass_imbalance", 1) >= 1e-6
            or record.get("draw_mass_fraction") != 0.065
            or record.get("feed_mass_fraction") != 0.00065
            or record.get("support_resistance_s_m") != 150666
            or record.get("inlet_offset_m") != 0.0045
            or record.get("inlet_diameter_m") != 0.0055
            or record.get("chamber_height_m") != 0.01925
            or not record.get("mesh_check_passed")
            or (record.get("geometry") == "CF042 toolkit geometry" and (
                record.get("channel_corner_radius_m") != 0.003
                or record.get("distributor_width_m") != 0.0295
            ))
        ):
            continue
        key = tuple(
            record.get(k) for k in ("flow_ml_min", "length_m", "width_m", "height_m")
        ) + (
            record.get("geometry", "CF042 reconstruction"),
            record.get("wall_normal_expansion", 5) if record.get("geometry") == "CF042 toolkit geometry" else None,
        )
        if key not in candidates or record.get("cells", 0) > candidates[key].get(
            "cells", 0
        ):
            candidates[key] = record
    original_ratios = sorted({r.get("wall_normal_expansion", 5) for r in candidates.values()
                             if r.get("geometry") == "CF042 toolkit geometry"})
    original_markers = {ratio: ("s", "^", "v", "P", "X")[i % 5]
                        for i, ratio in enumerate(original_ratios)}
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3), layout="constrained")
    defaults = dict(length=0.085, width=0.039, height=0.00225)
    for ax, dimension in zip(axes, ("length", "height", "width")):
        for flow, color in zip((5, 50, 500), (ORANGE, TEAL, INK)):
            reference = np.array(PAPER_2016_DIMENSIONS[dimension][flow])
            x = np.linspace(reference[:, 0].min(), reference[:, 0].max(), 120)
            values = []
            for value in x:
                p = {**defaults, dimension: value}
                kf = film_coefficient(flow, p["length"], p["width"], p["height"])
                values.append(
                    water_flux_mass_fractions(
                        0.065, 1.61111e-12, 8.33333e-8, 150666, 0.00065, kf
                    )
                    * 3.6e6
                )
            ax.plot(x * 1000, values, color=color, lw=1.7, label=f"{flow} mL/min")
            ax.scatter(
                reference[:, 0] * 1000, reference[:, 1], color=color, s=20, zorder=3
            )
            for record in candidates.values():
                if record["flow_ml_min"] != flow or any(
                    abs(record[other + "_m"] - defaults[other]) > 1e-10
                    for other in defaults
                    if other != dimension
                ):
                    continue
                ax.scatter(
                    record[dimension + "_m"] * 1000,
                    record["water_mass_flux_kg_m2_h"],
                    marker=(
                        original_markers[record.get("wall_normal_expansion", 5)]
                        if record.get("geometry") == "CF042 toolkit geometry"
                        else "D"
                    ),
                    facecolors="none",
                    edgecolors=color,
                    s=65,
                    lw=1.5,
                    zorder=4,
                )
        ax.axvline(defaults[dimension] * 1000, color=MUTED, lw=0.8, ls=":")
        ax.set(
            xlabel=f"{dimension.capitalize()} (mm)",
            ylim=(0, 17),
            title=f"Changing {dimension}",
        )
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Water mass flux (kg m$^{-2}$ h$^{-1}$)")
    axes[1].legend(frameon=False, ncol=3, loc="lower right", fontsize=8)
    symbols = [
        Line2D([], [], color=MUTED, label="Film model"),
        Line2D([], [], color=MUTED, marker="o", ls="", label="Archived CFD"),
    ]
    if any(r.get("geometry") != "CF042 toolkit geometry" for r in candidates.values()):
        symbols.append(
            Line2D(
                [],
                [],
                color=MUTED,
                marker="D",
                markerfacecolor="none",
                ls="",
                label="MembraneFoam",
            )
        )
    for ratio, marker in original_markers.items():
        symbols.append(Line2D([], [], color=MUTED, marker=marker, markerfacecolor="none",
                              ls="", label=f"Original toolkit, grading {ratio:g}"))
    axes[2].legend(handles=symbols, frameon=False, loc="lower right", fontsize=8)
    fig.suptitle(
        "2016 chamber dimensions • reference comparisons",
        fontsize=17,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.supxlabel(
        "Draw/feed salt mass fractions: 0.065 / 0.00065 • K = 150,666 s/m\n"
        "Lines: published film equations with reference density 1000 kg/m³. Dots: archived Figure 2 CFD values.",
        fontsize=9,
        color=MUTED,
    )
    save(fig, "2016-dimensions")


def chamber_results():
    """Select the finest validated 2012 grid at the published operating point."""
    selected = {}
    for source in RESULTS.glob("*.json"):
        record = json.loads(source.read_text())
        if not isinstance(record, dict):
            continue
        chamber = record.get("chamber")
        if chamber is None and source.stem in ("chamber-a", "chamber-b"):
            chamber = source.stem[-1].upper()  # Historical result exports.
        if (
            chamber not in ("A", "B")
            or record.get("geometry") != "2012 chamber"
            or record.get("flow_ml_min") != 50
            or record.get("draw_molar") != 1
            or not record.get("converged")
            or record.get("run_type") == "transient"
            or not record.get("mesh_check_passed")
            or record.get("relative_salt_imbalance", 1) >= 0.01
            or record.get("relative_mass_imbalance", 1) >= 1e-6
        ):
            continue
        if chamber not in selected or record.get("cells", 0) > selected[chamber][0].get(
            "cells", 0
        ):
            selected[chamber] = (record, source)
    return selected


def comparison():
    records = [(label, value[0]) for label, value in sorted(chamber_results().items())]
    if not records:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), layout="constrained")
    published = {
        "A": (5.46, 1.35, 5.64, 0.52, 1.44, 0.28),
        "B": (5.54, 1.37, 5.72, 0.40, 1.60, 0.39),
    }
    for ax, index, field, title, unit in zip(
        axes,
        (0, 1),
        ("water_mass_flux_kg_m2_h", "salt_flux_g_m2_h"),
        ("Water transport", "Reverse salt transport"),
        ("kg m$^{-2}$ h$^{-1}$", "g m$^{-2}$ h$^{-1}$"),
    ):
        for n, (label, r) in enumerate(records):
            row = published[label]
            ax.errorbar(
                n - 0.19,
                row[2 if index == 0 else 4],
                yerr=row[3 if index == 0 else 5],
                fmt="o",
                capsize=5,
                color=MUTED,
                label="2012 experiment" if n == 0 else None,
            )
            ax.plot(
                n,
                row[index],
                marker="D",
                ms=7,
                color=ORANGE,
                ls="",
                label="2012 published CFD" if n == 0 else None,
            )
            ax.plot(
                n + 0.19,
                r[field],
                marker="o",
                ms=9,
                color=TEAL,
                ls="",
                label="MembraneFoam" if n == 0 else None,
            )
            ax.annotate(
                f"{r[field]:.3f}",
                (n + 0.19, r[field]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=9,
                color=TEAL,
            )
        ax.set(
            xticks=range(len(records)),
            xticklabels=[f"Chamber {x[0]}" for x in records],
            ylabel=f"Flux ({unit})",
            title=title,
        )
        ax.set_xlim(-0.6, len(records) - 0.4)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    axes[0].legend(frameon=False, loc="lower left", fontsize=9)
    fig.suptitle(
        "2012 Table 2 • chamber validation",
        fontsize=17,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.supxlabel(
        "1 M NaCl draw • pure-water feed • 50 mL/min per compartment\nError bars are the published experimental variation; reconstructed geometry changes are documented.",
        fontsize=9,
        color=MUTED,
    )
    save(fig, "2012-flux-comparison")


def surface_maps():
    for label, (_, record_path) in sorted(chamber_results().items()):
        chamber = label.lower()
        source = record_path.with_name(record_path.stem + "-surface.csv")
        if not source.exists():
            continue
        data = np.genfromtxt(source, delimiter=",", names=True)
        draw = data[data["nz"] < -0.5]
        if not len(draw):
            continue
        # Draw velocity is scaled by density; convert to feed-equivalent mass flux.
        density = 997.1 + 694 * draw["mass_fraction"]
        flux = -draw["normal_velocity_m_s"] * density * 3600
        bulk = mass_fraction(1.0)
        feed = data[data["nz"] > 0.5]

        def coordinates(rows):
            return [tuple(row[k] for k in ("x_m", "y_m", "z_m")) for row in rows]

        paired = feed[paired_face_indices(coordinates(draw), coordinates(feed))]
        # 2012 Eq. 14: interface osmotic pressure from feed-side water velocity.
        interface = np.array(
            [
                row["normal_velocity_m_s"] / (1.222222222222e-12 * 80510000.0)
                + row["mass_fraction"]
                for row in paired
            ]
        )
        values = [draw["mass_fraction"] / bulk, interface / bulk, flux]
        fig, axes = plt.subplots(
            1, 3, figsize=(13.5, 6 if chamber == "a" else 4.8), layout="constrained"
        )
        for ax, val, title, cmap, unit in zip(
            axes,
            values,
            (
                "External polarization · Fig. 5",
                "Internal polarization · Fig. 6",
                "Water transport · Fig. 7",
            ),
            ("viridis", "cividis", "magma"),
            (
                "Surface / inlet mass fraction",
                "Interface / inlet mass fraction",
                "kg m$^{-2}$ h$^{-1}$",
            ),
        ):
            x = draw["x_m"] * 1000
            y = draw["y_m"] * 1000
            # Mirror the simulated half-width for a full physical surface view.
            artist = ax.scatter(
                np.r_[x, -x],
                np.r_[y, y],
                c=np.r_[val, val],
                s=7,
                marker="s",
                cmap=cmap,
                rasterized=True,
            )
            ax.set(
                aspect="equal",
                xlabel="Width coordinate (mm)",
                ylabel="Flow coordinate (mm)",
                title=title,
            )
            fig.colorbar(artist, ax=ax, shrink=0.75, label=unit)
        fig.suptitle(
            f"Chamber {chamber.upper()} • resolved membrane surface",
            fontsize=17,
            fontweight="bold",
            x=0.02,
            ha="left",
        )
        fig.supxlabel(
            "MembraneFoam • 2012 chamber geometry • half-width symmetry reflected for display.\n1 M draw / pure-water feed, 50 mL/min. Colors show computed values at membrane face centres.",
            fontsize=9,
            color=MUTED,
        )
        save(fig, f"2012-chamber-{chamber}-surface")


def channel_refinement():
    records = []
    for n in (40, 80, 160):
        source = RESULTS / f"channel-{n}.json"
        if source.exists():
            records.append((n, json.loads(source.read_text())))
    if len(records) < 2:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    cells = [2 * n * n for n, _ in records]
    fluxes = [r["water_mass_flux_kg_m2_h"] for _, r in records]
    axes[0].plot(cells, fluxes, "o-", color=TEAL, lw=2, ms=8)
    axes[0].set(
        xscale="log",
        xticks=cells,
        xticklabels=[f"{n:,}" for n in cells],
        xlabel="Cells",
        ylabel="Water mass flux (kg m$^{-2}$ h$^{-1}$)",
        title="Integral flux under refinement",
    )
    axes[0].minorticks_off()
    axes[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    for n, (_, record) in zip(cells, records):
        axes[0].annotate(
            f"{record['water_mass_flux_kg_m2_h']:.5f}",
            (n, record["water_mass_flux_kg_m2_h"]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )
    for (n, _), color in zip(records, (ORANGE, TEAL, INK)):
        source = RESULTS / f"channel-{n}-surface.csv"
        if not source.exists():
            continue
        data = np.genfromtxt(source, delimiter=",", names=True)
        draw = data[data["nz"] < -0.5]
        distance = (0.03 - draw["x_m"]) * 1000
        j = -draw["normal_velocity_m_s"] * (997.1 + 694 * draw["mass_fraction"]) * 3600
        order = np.argsort(distance)
        axes[1].plot(
            distance[order], j[order], lw=1.8, color=color, label=f"{2 * n * n:,} cells"
        )
    axes[1].set(
        xlabel="Distance from draw inlet (mm)",
        ylabel="Local water mass flux (kg m$^{-2}$ h$^{-1}$)",
        title="Resolved membrane profile",
    )
    axes[1].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    fig.suptitle(
        "Verification channel • spatial refinement",
        fontsize=17,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.supxlabel(
        "2012 membrane properties • 0.5 M draw / pure-water feed • 20 mL/min\nUniform-slot channel refinement; results apply to this geometry and operating point.",
        fontsize=9,
        color=MUTED,
    )
    save(fig, "channel-refinement")


def inlet_comparison():
    """Compare conserved, unspaced 3dBlock runs with the archived inlet studies."""
    from matplotlib.lines import Line2D
    from .paper_models import PAPER_2016_INLETS

    candidates = {}
    sources = {}
    for source in RESULTS.glob("*.json"):
        record = json.loads(source.read_text())
        if (
            not isinstance(record, dict)
            or record.get("geometry") != "3dBlock reconstruction"
        ):
            continue
        expected = dict(
            length_m=0.08,
            width_m=0.04,
            height_m=0.002,
            inlet_width_m=0.001,
            draw_mass_fraction=0.065,
            feed_mass_fraction=0.00065,
            support_resistance_s_m=150666,
        )
        if (
            not record.get("converged")
            or record.get("run_type") == "transient"
            or not record.get("mesh_check_passed")
            or record.get("spacer_count", 0)
            or record.get("relative_salt_imbalance", 1) >= 0.001
            or record.get("relative_mass_imbalance", 1) >= 1e-6
            or any(record.get(k) != v for k, v in expected.items())
        ):
            continue
        key = tuple(
            record.get(k) for k in ("flow_ml_min", "inlet_count", "inlet_angle_degrees")
        )
        if key[0] not in (5, 50, 500):
            continue
        if key not in candidates or record.get("cells", 0) > candidates[key].get(
            "cells", 0
        ):
            candidates[key] = record
            sources[key] = source
    selected = [
        r
        for r in candidates.values()
        if r.get("inlet_count") == 3 or r.get("inlet_angle_degrees") == 45
    ]
    if not selected:
        return  # A reference-only chart would not demonstrate reproduction.
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.5))
    fig.subplots_adjust(top=0.76, bottom=0.23, wspace=0.25)
    for ax, study, field, fixed, value, title in zip(
        axes,
        ("angle", "inlets"),
        ("inlet_angle_degrees", "inlet_count"),
        ("inlet_count", "inlet_angle_degrees"),
        (3, 45),
        ("Inlet angle · three inlets", "Inlet count · 45°"),
    ):
        for flow, color in zip((5, 50, 500), (ORANGE, TEAL, INK)):
            reference = np.asarray(PAPER_2016_INLETS[study][flow])
            ax.scatter(reference[:, 0], reference[:, 1], color=color, s=20)
            rows = [
                r for r in selected if r["flow_ml_min"] == flow and r[fixed] == value
            ]
            ax.scatter(
                [r[field] for r in rows],
                [r["water_mass_flux_kg_m2_h"] for r in rows],
                marker="D",
                s=60,
                facecolors="none",
                edgecolors=color,
                linewidths=1.5,
            )
        ax.set(
            title=title,
            xlabel="Angle (degrees)" if study == "angle" else "Number of inlets",
            ylabel="Water mass flux (kg m$^{-2}$ h$^{-1}$)",
        )
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        if study == "inlets":
            ax.set_xticks(range(1, 20, 2))
    handles = [
        Line2D([], [], color=c, marker="o", linestyle="none", label=f"{q} mL/min")
        for q, c in zip((5, 50, 500), (ORANGE, TEAL, INK))
    ]
    handles += [
        Line2D(
            [],
            [],
            color=MUTED,
            marker="D",
            markerfacecolor="none",
            linestyle="none",
            label="New CFD (finest accepted grid)",
        )
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        ncol=4,
        frameon=False,
    )
    fig.suptitle("2016 inlet geometry", fontsize=17, fontweight="bold")
    fig.supxlabel(
        "Filled circles: archived CFD · Open diamonds: new conserved solutions\n"
        "Agreement at individual points does not establish grid independence.",
        fontsize=9,
        color=MUTED,
    )
    save(fig, "2016-inlets")
    for key in sorted(candidates):
        if key[0] == 50 and key[2] == 45:
            inlet_surface(candidates[key], sources[key])


def spacer_comparison():
    """Compare source-layout count studies at their common base refinement."""
    from matplotlib.lines import Line2D
    from .paper_models import PAPER_2016_SPACERS

    expected = dict(length_m=0.08, width_m=0.04, height_m=0.002,
                    inlet_count=3, inlet_angle_degrees=45, inlet_width_m=0.001,
                    flow_ml_min=50, draw_mass_fraction=0.065,
                    feed_mass_fraction=0.00065, support_resistance_s_m=150666,
                    refinement=1)
    selected = {}
    selected_sources = {}
    baseline = None
    for source in sorted(RESULTS.glob("*.json")):
        record = json.loads(source.read_text())
        if not isinstance(record, dict):
            continue
        shape, placement = record.get("spacer_shape"), record.get("spacer_placement")
        count = record.get("spacer_count", 0)
        unspaced = (count == 0 and shape is None and placement is None
                    and record.get("spacer_placements") == []
                    and record.get("spacer_membrane_contact_count") == 0
                    and record.get("spacer_wall_contact_width_m") == 0)
        spaced = (shape in PAPER_2016_SPACERS
                  and placement in PAPER_2016_SPACERS[shape]
                  and count in range(2, 20, 2)
                  and record.get("spacer_layout") == "source"
                  and record.get("spacer_area_fraction") == 0.2)
        if (not (unspaced or spaced)
                or record.get("geometry") != "3dBlock reconstruction"
                or record.get("run_type") != "steady"
                or not record.get("converged")
                or not record.get("mesh_check_passed")
                or not record.get("relative_salt_imbalance", 1) < 0.001
                or not record.get("relative_mass_imbalance", 1) < 1e-6
                or any(record.get(k) != v for k, v in expected.items())):
            continue
        offset = record.get("spacer_profile_offset_m") or 0
        source_offset = (-1.5245078297943535e-5
                         if (shape, placement) == ("triangle", "middle") else 0)
        if not math.isclose(offset, source_offset, rel_tol=0, abs_tol=1e-12):
            continue
        schemes = re.sub(r"\s+", " ", record.get("configuration", {}).get(
            "system/fvSchemes", ""))
        if not all(entry in schemes for entry in (
                '"div(phi,U)" bounded Gauss upwind;',
                '"div(phi,m_A)" bounded Gauss upwind;',
                '"laplacian(mu,U)" Gauss cubic corrected;')):
            continue
        water = record.get("water_mass_flux_kg_m2_h", float("nan"))
        area = record.get("surface_integrals", {}).get("negative_z", {}).get("area_m2", 0)
        nominal = record.get("nominal_membrane_area_full_m2", 0)
        if not (math.isfinite(water) and water > 0 and area > 0 and nominal > 0):
            continue
        # Collected means use one symmetry half. Nominal area is full width.
        fraction = area * (2 if record.get("half_width") else 1) / nominal
        if not 0 < fraction <= 1 + 1e-9:
            continue
        values = (water, water * fraction, record.get("final_time", 0))
        if unspaced:
            if (math.isclose(fraction, 1, rel_tol=1e-9)
                    and math.isclose(nominal, expected["length_m"] * expected["width_m"],
                                     rel_tol=1e-9)
                    and (baseline is None or values[2] > baseline[2])):
                baseline = values
            continue
        key = (shape, placement, count)
        if key not in selected or record.get("final_time", 0) > selected[key][2]:
            selected[key] = values
            selected_sources[key] = (record, source)
    if not selected:
        return
    # N2 Changing and Down have identical generated geometry for each shape.
    for shape in PAPER_2016_SPACERS:
        if (shape, "down", 2) in selected:
            selected.setdefault((shape, "changing", 2), selected[(shape, "down", 2)])
        if baseline is not None:
            for placement in PAPER_2016_SPACERS[shape]:
                selected[(shape, placement, 0)] = baseline
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 7), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.79, bottom=0.19, left=0.07, right=0.99,
                        hspace=0.36, wspace=0.16)
    for row, shape in enumerate(("cylinder", "triangle")):
        for col, placement in enumerate(("middle", "up", "down", "changing")):
            ax = axes[row, col]
            reference = np.asarray(PAPER_2016_SPACERS[shape][placement])
            ax.scatter(reference[:, 0], reference[:, 1], color=MUTED, s=15)
            rows = sorted((n, values) for (s, p, n), values in selected.items()
                          if (s, p) == (shape, placement))
            ax.scatter([n for n, v in rows], [v[0] for n, v in rows],
                       marker="D", s=48, facecolors="none", edgecolors=TEAL)
            contact = [(n, v) for n, v in rows if abs(v[0] - v[1]) > 1e-7]
            ax.scatter([n for n, v in contact], [v[1] for n, v in contact],
                       marker="s", s=32, facecolors="none", edgecolors=ORANGE)
            ax.set(title=f"{shape.title()} · {placement.title()}",
                   xticks=(0, 6, 12, 18), xlim=(-1, 19))
            ax.grid(axis="y")
            ax.set_axisbelow(True)
            if row == 1:
                ax.set_xlabel("Spacers per channel")
            if col == 0:
                ax.set_ylabel("Water flux (kg m$^{-2}$ h$^{-1}$)")
    handles = [Line2D([], [], marker=m, color=c, markerfacecolor=fill,
                      linestyle="none", label=label)
               for m, c, fill, label in (
                   ("o", MUTED, MUTED, "Archived plotted means"),
                   ("D", TEAL, "none", "Current / active area"),
                   ("s", ORANGE, "none", "Current / nominal area (contact)"))]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.92),
               ncol=3, frameon=False)
    fig.suptitle("2016 spacer count · three inlets, 45°, 50 mL/min",
                 fontsize=17, fontweight="bold")
    fig.supxlabel("Source layout · equal-area fraction 0.2 · base refinement only\n"
                  "N0 is one shared unspaced baseline; N2 Changing aliases Down.\n"
                  "Only accepted current points shown. Historical normalization/settings and own-grid uncertainty remain unresolved.",
                  fontsize=9, color=MUTED)
    save(fig, "2016-spacers")
    middle = {n: selected_sources[("cylinder", "middle", n)]
              for n in (2, 18) if ("cylinder", "middle", n) in selected_sources}
    spacer_surface_maps(middle)


def spacer_surface_maps(records):
    """Compare accepted N2/N18 middle-cylinder fields on a shared flux scale."""
    if set(records) != {2, 18}:
        return
    fields = {}
    for count, (record, path) in records.items():
        source = path.with_name(path.stem + "-surface.csv")
        if not source.exists():
            return
        data = np.genfromtxt(source, delimiter=",", names=True)
        draw = data[data["nz"] < -0.5]
        if not len(draw):
            return
        material = record["configuration"]["constant/transportProperties"]
        coefficients = []
        for key in ("rho0", "rho_mACoeff"):
            entry = re.search(rf"\b{key}\s+(?:{key}\s+)?\[[^]]+\]\s+([^;]+);", material)
            if not entry:
                raise ValueError(f"Missing density coefficient {key} in {path}")
            coefficients.append(float(entry[1]))
        rho0, beta = coefficients
        flux = -draw["normal_velocity_m_s"] * rho0 * (1 + beta * draw["mass_fraction"]) * 3600
        area = np.sum(draw["area_m2"])
        if (not np.all(np.isfinite(flux)) or not np.all(draw["area_m2"] > 0)
                or not math.isclose(area, record["surface_integrals"]["negative_z"]["area_m2"],
                                    rel_tol=1e-9)
                or not math.isclose(np.sum(flux * draw["area_m2"]) / area,
                                    record["water_mass_flux_kg_m2_h"], rel_tol=1e-8)):
            raise ValueError(f"Surface export does not match accepted integral in {path}")
        x, y = draw["x_m"] * 1000, draw["y_m"] * 1000
        if record.get("half_width"):
            x, y, flux = np.r_[x, x], np.r_[y, -y], np.r_[flux, flux]
        fields[count] = (x, y, flux)
    low = min(values[2].min() for values in fields.values())
    high = max(values[2].max() for values in fields.values())
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), layout="constrained")
    for ax, count in zip(axes, (2, 18)):
        x, y, flux = fields[count]
        artist = ax.tripcolor(x, y, flux, shading="gouraud", cmap="magma",
                              vmin=low, vmax=high, rasterized=True)
        mean = records[count][0]["water_mass_flux_kg_m2_h"]
        ax.set(aspect="equal", xlim=(0, 80), ylim=(-20, 20),
               xlabel="Length (mm)", ylabel="Width (mm)",
               title=f"{count} middle cylinders · mean {mean:.6f} kg/(m² h)")
    fig.colorbar(artist, ax=axes, shrink=0.85,
                 label="Local water mass flux (kg m$^{-2}$ h$^{-1}$)")
    fig.suptitle("2016 middle-cylinder spacers · membrane water flux",
                 fontsize=17, fontweight="bold")
    fig.supxlabel("Three inlets · 45° · 50 mL/min · source layout · area fraction 0.2 · base refinement\n"
                  "Accepted face-centre values interpolated; half-width reflected. Historical fields and grid uncertainty remain unverified.",
                  fontsize=9, color=MUTED)
    save(fig, "2016-spacer-surfaces")


def inlet_surface(record, record_path):
    """Render an accepted inlet-count case; reflect its half-width symmetry."""
    source = record_path.with_name(record_path.stem + "-surface.csv")
    if not source.exists():
        return
    data = np.genfromtxt(source, delimiter=",", names=True)
    draw, feed = data[data["nz"] < -0.5], data[data["nz"] > 0.5]
    if not len(draw) or not len(feed):
        return
    flux = -draw["normal_velocity_m_s"] * (997.1 + 694 * draw["mass_fraction"]) * 3600
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), layout="constrained")
    for ax, rows, values, title, unit, cmap in zip(
        axes,
        (draw, feed, draw),
        (draw["mass_fraction"] / record["draw_mass_fraction"],
         feed["mass_fraction"] / record["feed_mass_fraction"], flux),
        ("Draw-side concentration", "Feed-side concentration", "Local water flux"),
        ("Surface / inlet mass fraction", "Surface / inlet mass fraction", "kg m$^{-2}$ h$^{-1}$"),
        ("viridis", "cividis", "magma"),
    ):
        x, y = rows["x_m"] * 1000, rows["y_m"] * 1000
        artist = ax.tripcolor(np.r_[x, x], np.r_[y, -y], np.r_[values, values],
                              shading="gouraud", cmap=cmap, rasterized=True)
        ax.set(aspect="equal", xlim=(0, 80), ylim=(-20, 20),
               xlabel="Length (mm)", ylabel="Width (mm)", title=title)
        fig.colorbar(artist, ax=ax, shrink=0.8, label=unit)
    count = record["inlet_count"]
    fig.suptitle(f"2016 chamber with {count} inlets · membrane fields", fontsize=17, fontweight="bold")
    fig.supxlabel(f"{count} inlets · 45° · 50 mL/min · reconstructed geometry\n"
                  "Interpolated face-centre values; half-width reflected. Grid independence remains unverified.",
                  fontsize=9, color=MUTED)
    save(fig, "2016-inlet-surface" if count == 3 else f"2016-inlet-{count}-surface")


def main():
    global RESULTS, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=RESULTS,
        help="Collected case data (default: runs/results)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUT,
        help="Figure directory (default: docs/figures)",
    )
    args = parser.parse_args()
    RESULTS, OUT = args.results, args.output
    if not RESULTS.exists():
        print(f"No CFD records at {RESULTS}; generating analytical figures only.")
    if all((RESULTS / f"channel-{n}.json").exists() for n in (40, 80, 160)):
        from .results import refinement

        refinement(RESULTS)
    analytical()
    inlet_comparison()
    spacer_comparison()
    comparison()
    surface_maps()
    channel_refinement()
