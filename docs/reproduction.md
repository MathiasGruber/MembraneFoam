# Paper comparisons

## 2012: experimental validation and spatial polarization

Source: Gruber et al., DOI 10.3390/membranes2040764. Parameters in section 3.1:

| Parameter | Published | SI value used |
|---|---:|---:|
| A | 0.44 +/- 0.05 LMH/bar | 1.222222222222e-12 m/(Pa s) |
| B | 0.087 +/- 0.018 LMH | 2.416666666667e-8 m/s |
| K | 0.72 +/- 0.23 s/um | 720000 s/m |

Table 2 and Figures 5-7 use 1 M draw, pure-water feed, and 50 mL/min per
compartment. Simulations use half-width symmetry, so each simulated compartment
receives 25 mL/min. Molarity is converted with the concentration-dependent density.

| Chamber | Published CFD water, kg/(m2 h) | Published CFD salt, g/(m2 h) | Experimental water | Experimental salt |
|---|---:|---:|---:|---:|
| A | 5.46 | 1.35 | 5.64 +/- 0.52 | 1.44 +/- 0.28 |
| B | 5.54 | 1.37 | 5.72 +/- 0.40 | 1.60 +/- 0.39 |

The chamber generators derive from the September 2012 `createSMTCdict.py` (A)
and `createFO21Dict.py` (B). Dimensions match the published chambers; exact
identity with the publication meshes is unverified.

Chamber A uses millimetre scaling (`convertToMeters=0.001`), omits four unused
curved edges, and sets the inlet connection offset to 0.05 mm instead of
0.025 mm to avoid sliver faces.

The generators preserve the original geometry construction. Mesh assembly uses
`mirrorMesh -overwrite`, face-zone selection, and `createBaffles`. The inlet/outlet
patches are selected geometrically rather than split using hard-coded face counts.
`checkMesh` must pass before assigning fields and running a case.

## 2012 simulation results

Both meshes pass `checkMesh`. Eight-rank CPU results:

| Chamber | Cells (simulated half-width) | Water, kg/(m² h) | Salt, g/(m² h) | Salt imbalance / membrane transfer |
|---|---:|---:|---:|---:|
| A | 326,720 | 5.833361 | 1.436791 | 0.036% |
| B | 553,560 | 5.615712 | 1.383162 | 0.150% |

Water flux differs from published CFD by +6.8% (A) and +1.4% (B); salt flux differs
by +6.4% and +1.0%. All four predictions fall within the published experimental
ranges. Geometry and discretization differences prevent a direct comparison
of the historical numerical fields.

The cases use linear-upwind momentum convection, van Leer salt convection,
two non-orthogonal pressure corrections, salt relaxation 0.2, and linear pressure
tolerance 1e-12. SIMPLE residual limits are 1e-6 for pressure/velocity and 1e-9
for salt, with the additional integral salt-conservation stopping condition.
Pressure/velocity relaxation is 0.15/0.3 for A and 0.3/0.7 for B.
Chamber A retains 100 cells at the coarsest pressure level and solves that level
directly. It uses DIC–Gauss–Seidel smoothing for flow initialization, then restores
the steady solver settings.
The example runner exports dictionaries, residuals, and membrane surface values
to the chosen results directory.

Small interior mass-fraction undershoots remain: minima are approximately
-3.7e-8 (A) and -3.0e-8 (B), compared with a draw mass fraction of 0.0564.
They are recorded without clipping. Membrane concentrations remain admissible.
Three-dimensional chamber grid independence is unverified. Surface maps have
not been quantitatively compared with the published spatial distributions.

Use `--refinement` to increase the block cell counts while preserving dimensions
and grading:

```bash
for refinement in 1 1.25 1.5; do
    for chamber in a b; do
        case="chamber-$chamber-refinement-$refinement"
        python3 examples/run.py "chamber-$chamber" "runs/$case" \
            --refinement "$refinement" --ranks 8 --results "runs/results/$case"
    done
done
```

At refined port intersections, the runner checks the quarter mesh and uses
OpenFOAM's `collapseEdges` to remove identified sliver faces before mirroring.
The edge-length threshold is 10 nm; the collapse filter's minimum volume is
1e-18 m³ to accommodate the millimetre-scale channels. The assembled mesh must
still pass the standard `checkMesh` checks. Refinement factors can interact with
the legacy port intersections, so a finer mesh is not automatically a valid mesh.
Passing mesh checks alone does not establish flux or surface-field convergence.
The 2012 comparison and surface plots select the finest converged, conserved
record for each chamber at the published operating point.

Figures 5–7 are recreated as surface maps. External polarization is membrane
surface mass fraction divided by inlet draw mass fraction. Internal polarization
uses the support/active-layer interface concentration from Eq. 14, divided by
the same inlet value. Water maps show mass flux. Colour scales are independent
between chambers; consult the labelled ranges when comparing them. The plots
reflect the simulated half-width and display face-centre values without interpolation.

## 2016 chamber dimensions

Source: Gruber, Aslak & Hélix-Nielsen, DOI 10.1016/j.seppur.2015.12.017.
Figure 2 compares 5, 50 and 500 mL/min with the rectangular-channel film model,
Eqs. 4 and 9–11. Reference CFD values are transcribed from the authors' Figure 2
numerical tables in `paper_models.py`; they are benchmark inputs.

| Parameter | Value |
|---|---:|
| Membrane length / width | 85 / 39 mm |
| Channel height per compartment | 2.25 mm |
| Total compartment height | 19.25 mm |
| Inlet diameter / edge offset | 5.5 / 4.5 mm |
| Draw / feed salt mass fraction | 0.065 / 0.00065 kg/kg |
| Water permeability A | 1.61111e-12 m/(Pa s) |
| Salt permeability B | 8.33333e-8 m/s |
| Support resistance K | 150666 s/m |

The original `Jw_vs_dimensions/createPlots.py` specifies these membrane
parameters and concentrations. The concentrations also appear explicitly in the
[historical case initialization](../tests/simple3dMesh/system/setFieldsDict).
The analytical calculation uses constant diffusivity 1.45e-9 m²/s and reference
density 1000 kg/m³ for conversion to mass flux. Independent samples from the original
analytical curves agree within 0.06 kg/(m² h), the figure's digitization tolerance.
The CFD solvers retain concentration-dependent density and diffusivity.

Table 1 prints K=6.64 s/µm, which converts to 6,640,000 s/m. This conflicts with
the original plotting code and case dictionaries. It is not used for the
reference comparison: even without external concentration polarization it
limits flux to about 1.46 kg/(m² h), below the archived CFD range of 8.80–13.22.
The recovered value is 0.150666 s/µm; its reciprocal is 6.6372 µm/s, which
rounds to the table's 6.64 with inverse units. This is consistent with a
reciprocal/unit error in the table, although the original numerical inputs
do not establish how the error arose. The reference cases use those inputs
rather than fitting K to the published fluxes.

The CF042 generator uses half-width symmetry, paired membrane faces, rounded
channels and distributors, and circular ports. Each compartment receives half
the specified whole-module flow. Feed and draw enter opposite ends. The
reconstruction assumes a channel corner radius of half the inlet diameter,
distributor half-span `width/2 - diameter/2`, and distributor top height
`chamber_height - diameter/2`; a curved block transition connects it to the pipe.
These details are not uniquely specified by the published drawings. This is a
reconstruction, and exact identity with the original mesh is unverified.

Default tangential spacing is 0.8 mm with 24 cells through each channel and
wall-normal expansion ratio 100, giving approximately 4 µm at the membrane.
Use `--spacing 0.0004 --layers 48` for a finer grid. Every mesh must pass
`checkMesh`. Solver stopping requires salt residual below 1e-9 and salt
imbalance below 0.1% of membrane transfer. The plotter includes new CFD points
only after convergence and conservation checks; original CFD points are always
labelled separately.

Converged CF042 results:

| Flow (mL/min) | Cells | Water, kg/(m² h) | Salt, g/(m² h) | Difference from archived water flux | Salt imbalance / transfer |
|---:|---:|---:|---:|---:|---:|
| 5 | 158,912 | 8.382253 | 5.392870 | −4.7% | 0.0114% |
| 50 | 158,912 | 11.180289 | 7.196734 | −1.7% | 0.0046% |
| 50 | 396,340 | 11.184071 | 7.199158 | −1.6% | <0.1000% |

Relative mass imbalance is below 4e-12. The two 50 mL/min grids differ by 0.034%
in water flux and 0.081% in the draw-surface concentration integral. A third
grid is needed to establish chamber grid independence; the difference from the
archived CFD results remains unresolved.

Initialize a refined mesh from an existing solution to reduce startup cost:

```bash
python3 examples/run.py cf042 runs/cf042-50-fine --ranks 8 \
    --spacing 0.0004 --layers 48 --initialize-from runs/cf042-50 \
    --results runs/results/cf042-50-fine
```

The source must have the same geometry and operating conditions. The runner uses
`mapFields` with the latest source time, supports decomposed sources, and resets
the coincident membrane boundary guesses before solving. This initialization
does not change the convergence or conservation requirements.

## 2016 inlet geometry

The `3dblock` example generates the 80 × 40 mm chamber with a 2 mm channel
on each side of the membrane. Its rounded ends and 1 mm rectangular ports
follow the [archived quarter mesh](../tests/simple3dMesh/constant/polyMesh/blockMeshDict).
`--inlets` selects an odd whole-width count from 1 to 19; `--angle` selects
0–90° relative to the membrane. The default is three inlets at 45° and
50 mL/min. Feed and draw flow in opposite directions, with half-width symmetry.
The membrane parameters and concentrations are those listed above.

```bash
python3 examples/run.py 3dblock runs/block-3 --inlets 3 --angle 45 --ranks 8
python3 examples/run.py 3dblock runs/block-3-fine --inlets 3 --angle 45 \
    --refinement 1.25 --ranks 8
```

The port axes rotate about a point 0.5 mm toward the inlet from the rounded-end centre
and 0.5 mm above the membrane. Port faces meet the curved chamber wall directly.
The mesh partitions that wall at both port lips, avoiding thin cut cells.
`--refinement` increases cell counts in every direction while retaining the
geometry and wall-normal grading. Each generated mesh must pass `checkMesh`.
A passing mesh establishes numerical quality; agreement with the published
inlet-count and angle curves requires converged, mesh-resolved solutions.

Generate inlet-angle and inlet-count comparisons with the same runner:

```bash
for angle in 0 45 90; do
    python3 examples/run.py 3dblock "runs/inlet-angle-$angle" --angle "$angle" \
        --ranks 8 --results "runs/results/inlet-angle-$angle"
done
for count in 1 7; do
    python3 examples/run.py 3dblock "runs/inlet-count-$count" --inlets "$count" \
        --ranks 8 --results "runs/results/inlet-count-$count"
done
python3 examples/plot.py --results runs/results
```

The plotter writes `2016-inlets.svg` when eligible inlet-study results exist.
Filled circles are the authors' numerical reference tables; open diamonds are
new solutions. It selects the finest accepted mesh at each operating point,
excludes spacer and transient cases, and requires salt imbalance below 0.1%
of membrane transfer. Use `--flow` to repeat the studies at 5 and 500 mL/min.

At 50 mL/min, the accepted three-inlet solutions are:

| Angle | Cells | Water flux, kg/(m² h) | Archived flux | Salt imbalance |
| --- | ---: | ---: | ---: | ---: |
| 30° | 1,066,000 | 12.018934 | 11.90 | 0.09996% |
| 45° | 1,066,000 | 11.897860 | 11.90 | 0.09814% |
| 50° | 1,066,000 | 11.861254 | 11.90 | 0.09912% |
| 60° | 1,066,000 | 11.823157 | 11.89 | 0.09790% |
| 70° | 1,030,800 | 11.817270 | 11.89 | 0.03199% |
| 80° | 1,014,800 | 11.793435 | 11.89 | 0.08206% |
| 90° | 1,014,800 | 11.758823 | 11.89 | 0.04556% |

Salt imbalance is relative to membrane transfer; relative mass imbalance is
below 5e-12 for every listed case. The accepted 45°, 50 mL/min inlet-count points are:

| Inlet count | Cells | Water flux, kg/(m² h) | Archived flux | Salt imbalance |
| --- | ---: | ---: | ---: | ---: |
| 5 | 1,078,476 | 11.810781 | 11.91 | 0.00176% |
| 7 | 1,080,400 | 11.776103 | 11.92 | 0.00945% |
| 9 | 1,087,600 | 11.765107 | 11.92 | 0.09999% |
| 11 | 1,084,248 | 11.755738 | 11.92 | 0.04414% |
| 15 | 1,109,200 | 11.739117 | 11.92 | 0.07290% |

Their relative mass imbalances are 1.5e-12, 3.0e-12, 5.6e-12, 2.5e-12, and
1.7e-12. The fluxes lie 0.8% to 1.5% below the archived values, with the offset
growing from five to fifteen inlets.
The 10° solution reached the time limit without meeting the 0.1% salt-conservation gate: the relative salt imbalance decayed to 0.115% of the membrane transfer by about a fifth of the run and then plateaued for the remainder, ending at 0.1145%. At the limit the water flux is 12.197 kg/(m² h), 2.4% above the archived 11.91, and the relative mass imbalance is below 3e-13. It is therefore excluded from the accepted table, and the archived 10° point is not reproduced. The 0° solution is numerically unstable: every steady-state trial terminates during the early transient with nonphysical values (negative salt mass fraction on the membrane), and low-order discretization variants remain far above the 0.1% salt-conservation gate when their iteration budget ends. The archived 0° point is therefore not reproduced.
These points do not establish grid independence or reproduce
the complete inlet studies. At 50 mL/min and 45°, the plotter also writes
membrane-field maps for each accepted inlet count.
The 3dBlock runner uses a linear pressure tolerance of 1e-11 and relative
tolerance of 0.001. This avoids driving GAMG below the attainable residual on
the graded inlet meshes. Residual convergence and the 0.1% salt-conservation
requirement must both pass before a result is accepted.

### Middle-cylinder spacers

`--spacers` adds an even number of transverse cylinders, uniformly spaced along
each channel. Their centres lie halfway between the membrane and the opposite
wall, at longitudinal positions `length * i / (count + 1)`. The cylinders span
the full width. `--spacer-radius` is an explicit reconstruction parameter: the
archived cross-section drawing does not specify a numerical radius.

```bash
python3 examples/run.py 3dblock runs/block-2-spacers --spacers 2 \
    --spacer-radius 0.0005 --ranks 8
```

For this example radius of 0.5 mm, each cylinder displaces
`π * radius² * width` of fluid per compartment. The mesh grades toward both
the cylinder and channel walls. Use `--refinement` to test spatial convergence
and vary the radius to assess geometric sensitivity. These are reconstructed
middle-cylinder cases; the triangle and wall-contact configurations are not
implemented, and agreement with the published spacer curves is unverified.

The surface CSV includes `area_m2`. Run summaries report the area, concentration
integral, and area-weighted mean under `surface_integrals`, grouped by the sign
of the face's z-normal. In these FO examples, `negative_z` is the draw side.
Use its concentration integral for chamber grid comparisons, as in the paper;
an unweighted mean of face values is inappropriate on the graded meshes.

## Transient integration

Use `--transient` and `--time-step` to run PISO. SIMPLE relaxation factors are
cleared, and `--outer-correctors` selects the coupling iterations per step
(default: four). `--non-orthogonal-correctors` controls repeated pressure solves
for mesh non-orthogonality (default: two). `--momentum-relaxation` optionally
damps the momentum equation with a factor in `(0, 1]` (default: one). Pressure
and solute equations remain unrelaxed. Check sensitivity to the time step,
corrector counts, and relaxation; the settings must preserve conservation and
the resolved surface fields:

```bash
python3 examples/run.py cf042 runs/cf042-500-transient --flow 500 --ranks 8 \
    --transient 60 --time-step 0.00025 --initialize-from runs/cf042-500
```

The 2016 paper followed steady solutions with at least 60 seconds of transient
flow to investigate time dependence. The initial state is mapped onto a new case
and its physical clock starts at zero. The duration must be a multiple of the
fixed time step. Choose it using the reported Courant numbers, and verify
time-step sensitivity. The runner writes approximately twenty checkpoints,
always saves the final state, and exports
`result-history.csv` (or the chosen result stem), containing sampled salt and mass
balances including accumulation. The FO flux at each sample uses the final PISO
update preceding that audit.

Collection requires every sampled salt-balance error below 0.1% of membrane
transfer and mass-balance error below 1e-6 of throughput. Concentration bounds
must remain physical at every completed time step, including steps between
saved checkpoints. These checks establish
conservation, not stationarity. Transient summaries retain `converged: false` and
are excluded from steady paper comparisons.
Steady runs that stop at their time limit without reporting SIMPLE convergence are likewise rejected by the collector and are not reported as reproduced steady results. Establish a suitable averaging
window, temporal convergence, and spatial convergence before comparing an
unsteady simulation with a published steady result.

## Independent channel refinement

The uniform-slot verification channel uses 0.5 M draw, pure-water feed and
20 mL/min. Both resolved directions are refined by two at each level; the width
is an empty direction. These results apply to the uniform channel.

| Cells | Water flux, kg/(m2 h) | Relative mass imbalance |
|---:|---:|---:|
| 3,200 | 4.11873389 | 2.42e-14 |
| 12,800 | 4.11788571 | 1.12e-14 |
| 51,200 | 4.11747221 | 5.85e-14 |

The observed order is approximately 1.04. `examples/plot.py`
computes the fine-grid GCI with refinement ratio 2 and safety factor 1.25;
the result is approximately 0.012%. This estimates discretization uncertainty
for this channel and operating point, not model uncertainty or 3D chamber GCI.
The two smaller cases ran serially; the fine case used four CPU ranks. The
separate serial/MPI regression checks agreement at the same mesh resolution.

The GCI calculation follows the three-grid procedure in the
[NASA spatial-convergence tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).
It is an estimate for the reported numerical tolerances, not a guarantee of error
bounds. The medium-to-fine flux change is 0.010%; the last 100 logged flux
evaluations span 0.000012% on the fine mesh.

## Generate the figures

Run from the repository root with OpenFOAM and MembraneFoam loaded. Each case
requires a fresh directory. Choose the MPI rank count for the available CPU resources.

```bash
for chamber in a b; do
    python3 examples/run.py "chamber-$chamber" "runs/chamber-$chamber" --ranks 8 \
        --results "runs/results/chamber-$chamber"
done
for n in 40 80 160; do
    python3 examples/run.py channel "runs/channel-$n" --nx "$n" --nz "$n" \
        --results "runs/results/channel-$n"
done
python3 -m pip install -r examples/requirements.txt
python3 examples/plot.py --results runs/results --output docs/figures
```

The plotter creates the 2012 comparison and surface maps from chamber records,
the refinement figure from channel records, and the 2016 comparison from
the film equations, archived reference values and eligible CF042 records.
Missing CFD records are skipped. SVG outputs are saved
to `--output`; collected JSON/CSV data and full cases remain in ignored `runs/`.
