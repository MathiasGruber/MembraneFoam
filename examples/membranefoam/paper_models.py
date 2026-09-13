"""Independent scalar reference calculations for the published analytical models.

These evaluate equations, not CFD. Concentrations at the membrane may differ
from inlet concentrations; the optional film model approximates that difference.
"""

import math


def mass_fraction(molar):
    concentration = molar * 58.44
    return 2 * concentration / (997.1 + math.sqrt(997.1**2 + 4 * 694 * concentration))


def water_flux(draw_molar, a, b, k, feed_molar=0.0, film_coefficient=None):
    return water_flux_mass_fractions(
        mass_fraction(draw_molar), a, b, k, mass_fraction(feed_molar), film_coefficient
    )


def water_flux_mass_fractions(draw, a, b, k, feed=0.0, film_coefficient=None):
    """Water velocity in m/s for bulk salt mass fractions (kg/kg)."""
    if not (
        all(math.isfinite(v) for v in (draw, feed, a, b, k))
        and 0 <= feed <= draw <= 0.09
        and a > 0
        and b >= 0
        and k >= 0
    ):
        raise ValueError("Invalid membrane parameters or bulk mass fractions")
    if film_coefficient is not None and (
        not math.isfinite(film_coefficient) or film_coefficient <= 0
    ):
        raise ValueError("Film coefficient must be positive")
    if draw == feed:
        return 0.0
    phi = 80510000.0
    lo, hi = 0.0, a * phi * (draw - feed)
    if hi < 0:
        raise ValueError("AL-FS reference requires draw >= feed")

    # Independent log-form residual for positive K and B.
    def residual(j):
        e = j / film_coefficient if film_coefficient else 0.0
        pm_feed = phi * feed * math.exp(min(e, 700))
        pm_draw = phi * draw * math.exp(-e)
        if k == 0:
            return j - a * (pm_draw - pm_feed)
        return k * j - math.log((b + a * pm_draw) / (b + j + a * pm_feed))

    for _ in range(100):
        mid = (lo + hi) / 2
        if residual(mid) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def film_coefficient(flow_ml_min, length, width, height, diffusivity=1.45e-9):
    if not all(
        math.isfinite(v) and v > 0
        for v in (flow_ml_min, length, width, height, diffusivity)
    ):
        raise ValueError("Flow, dimensions and diffusivity must be finite and positive")
    speed = flow_ml_min * 1e-6 / 60 / (width * height)
    return 1.62 * (speed * diffusivity**2 / (2 * height * length)) ** (1 / 3)


# Reference CFD values accompanying Figure 2, Gruber et al. (2016).
# Coordinates are metres; water fluxes are kg/(m² h).
PAPER_2016_DIMENSIONS = {
    "length": {
        5: [
            [0.05, 9.035],
            [0.085, 8.8],
            [0.1, 8.6765],
            [0.2, 7.9872],
            [0.3, 7.5161],
            [0.4, 7.1699],
            [0.5, 6.91525],
        ],
        50: [
            [0.05, 11.54],
            [0.085, 11.37],
            [0.1, 11.28],
            [0.2, 10.74],
            [0.3, 10.34],
            [0.4, 10.04],
            [0.5, 9.8],
        ],
        500: [
            [0.05, 13.33],
            [0.085, 13.22],
            [0.1, 13.15],
            [0.2, 12.79],
            [0.3, 12.5239],
            [0.4, 12.3104],
            [0.5, 12.1348],
        ],
    },
    "height": {
        5: [
            [0.004, 7.58414],
            [0.0035, 7.8737],
            [0.002875, 8.2943],
            [0.00225, 8.8],
            [0.001625, 9.42957],
            [0.001, 10.26],
        ],
        50: [
            [0.004, 10.369],
            [0.0035, 10.61],
            [0.002875, 10.965],
            [0.00225, 11.37],
            [0.001625, 11.88],
            [0.001, 12.54],
        ],
        500: [
            [0.004, 12.5741],
            [0.0035, 12.7369],
            [0.002875, 12.9623],
            [0.00225, 13.22],
            [0.001625, 13.53],
            [0.001, 13.925],
        ],
    },
    "width": {
        5: [
            [0.039, 8.8],
            [0.05, 8.58],
            [0.06, 8.37],
            [0.07, 8.17],
            [0.08, 7.98],
            [0.09, 7.79],
            [0.1, 7.62],
            [0.11, 7.45],
            [0.12, 7.27],
            [0.13, 7.12],
            [0.14, 6.95],
            [0.15, 6.8],
        ],
        50: [
            [0.039, 11.37],
            [0.05, 11.21],
            [0.06, 11.04],
            [0.07, 10.89],
            [0.08, 10.74],
            [0.09, 10.58],
            [0.1, 10.43],
            [0.11, 10.3],
            [0.12, 10.16],
            [0.13, 10.02],
            [0.14, 9.88],
            [0.15, 9.73],
        ],
        500: [
            [0.039, 13.22],
            [0.05, 13.11],
            [0.06, 12.99],
            [0.07, 12.88],
            [0.08, 12.78],
            [0.09, 12.67],
            [0.1, 12.58],
            [0.11, 12.47],
            [0.12, 12.38],
            [0.13, 12.28],
            [0.14, 12.17],
            [0.15, 12.07],
        ],
    },
}


# Authors' numerical tables: CFD_ArticleData/Jw_vs_{angle,inlets}/data.
# Angles in degrees; water mass flux in kg/(m² h).
PAPER_2016_INLETS = {
    "angle": {
        5: [
            [0.0, 9.48],
            [10.0, 9.46],
            [20.0, 9.47],
            [30.0, 9.47],
            [40.0, 9.47],
            [50.0, 9.46],
            [60.0, 9.46],
            [70.0, 9.46],
            [80.0, 9.46],
            [90.0, 9.46],
        ],
        50: [
            [0.0, 11.91],
            [10.0, 11.91],
            [20.0, 11.91],
            [30.0, 11.9],
            [40.0, 11.9],
            [50.0, 11.9],
            [60.0, 11.89],
            [70.0, 11.89],
            [80.0, 11.89],
            [90.0, 11.89],
        ],
        500: [
            [0.0, 13.58],
            [10.0, 13.57],
            [20.0, 13.57],
            [30.0, 13.56],
            [40.0, 13.55],
            [50.0, 13.55],
            [60.0, 13.54],
            [70.0, 13.54],
            [80.0, 13.54],
            [90.0, 13.54],
        ],
    },
    "inlets": {
        5: [
            [1.0, 9.38],
            [3.0, 9.46],
            [5.0, 9.48],
            [7.0, 9.5],
            [9.0, 9.51],
            [11.0, 9.51],
            [13.0, 9.5],
            [15.0, 9.5],
            [17.0, 9.49],
            [19.0, 9.48],
        ],
        50: [
            [1.0, 11.85],
            [3.0, 11.9],
            [5.0, 11.91],
            [7.0, 11.92],
            [9.0, 11.92],
            [11.0, 11.92],
            [13.0, 11.92],
            [15.0, 11.92],
            [17.0, 11.92],
            [19.0, 11.91],
        ],
        500: [
            [1.0, 13.52],
            [3.0, 13.55],
            [5.0, 13.55],
            [7.0, 13.55],
            [9.0, 13.56],
            [11.0, 13.56],
            [13.0, 13.56],
            [15.0, 13.55],
            [17.0, 13.55],
            [19.0, 13.54],
        ],
    },
}
