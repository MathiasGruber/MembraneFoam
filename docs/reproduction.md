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
| B | 553,560 | 5.615727 | 1.383166 | 0.08230% |
| B | 1,140,744 | 5.627176 | 1.385987 | 0.0971% |
| B | 1,926,360 | 5.628058 | 1.386204 | 0.09989% |

On the finest accepted meshes, water flux differs from published CFD by +6.8%
(A) and +1.6% (B); salt flux differs by +6.4% and +1.2%. All four predictions
fall within the published experimental ranges. Geometry and discretization
differences prevent a direct comparison of the historical numerical fields.

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
Chamber B has three accepted meshes under the same 0.1% integral
salt-imbalance gate. The strict coarse run replaces the earlier result obtained
under a 1% gate. Water flux changes +0.204% from the coarse to the 1.25× mesh,
then +0.01567% on the 1.5× mesh. The draw-surface concentration integrals are:

| B cells | Draw-side concentration integral, m² |
|---:|---:|
| 553,560 | 1.7574121e-5 |
| 1,140,744 | 1.7645680e-5 |
| 1,926,360 | 1.7651813e-5 |

Independent CSV integration confirms both membrane areas and integrals on all
three meshes. The two finest integrals differ by +0.03476%. Using effective
spacing proportional to `cells**(-1/3)` and the nonuniform refinement ratios,
the apparent orders are about 10.2 for water flux and 9.8 for the draw integral.
These high fitted orders do not establish an asymptotic error estimate; a
fourth mesh and iterative-error sensitivity are still needed before treating
a Richardson extrapolation or GCI as validated. The comparison and maps use
the finest accepted mesh. Chamber A grid independence and quantitative agreement
with the published spatial distributions remain unverified.

The [2012 paper, §2.3, Eq. (7)](https://doi.org/10.3390/membranes2040764), prints
`GCI_coarse = 3*abs(e)*R**3/(R**3 - 1)` and defines `R` as the fine/coarse
cell-count ratio. The conventional GCI uses a grid-spacing ratio and a
convergence order; in three dimensions the effective spacing ratio is
`(N_fine/N_coarse)**(1/3)`, as described in
[NASA's spatial-convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).
The printed definition differs from this prescription; its historical
implementation is unverified. Applying the printed cell-count formula with an
explicitly coarse-normalized concentration-integral difference gives 1.36%
for 553,560 versus 1,926,360 cells, or 0.132% for the two finest meshes. These
are comparisons with the printed criterion, not validated uncertainty bounds.

The finest B study uses PCG with a GAMG preconditioner after standalone GAMG
pressure solves stalled near their 1e-12 absolute tolerance. This changes the
linear solution method, preserving the equations, tolerances and acceptance
gates. OpenFOAM v2606 requires the preconditioner's controls in a nested
`preconditioner` dictionary. To reproduce this configuration from a fresh case:

```bash
case=runs/chamber-b-refinement-1.5
python3 examples/run.py chamber-b "$case" --refinement 1.5 --end 60000 \
    --salt-imbalance 0.001 --setup-only
foamDictionary "$case/system/fvSolution" -entry solvers.p -set '{
    solver PCG;
    preconditioner {
        preconditioner GAMG;
        smoother DICGaussSeidel;
        tolerance 1e-6;
        relTol 0.1;
        nVcycles 2;
    }
    tolerance 1e-12;
    relTol 0.001;
    maxIter 2000;
}'
cat > "$case/system/decomposeParDict" <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 8;
method simple;
simpleCoeffs { n (8 1 1); delta 0.001; }
constraints { membranePairs { type preserveBaffles; } }
EOF
decomposePar -case "$case" > "$case/log.decomposePar" 2>&1
mpirun -np 8 simpleSaltTransport -parallel -case "$case" \
    > "$case/log.solver" 2>&1
python3 examples/run.py chamber-b "$case" --collect-only \
    --results runs/results/chamber-b-refinement-1.5
```

Use `--salt-imbalance 0.001` to require the stricter integral gate during
steady solving and collection. The requested limit is stored in the result
metadata and cannot be replaced by the collector’s default 1% limit.

Use `--refinement` to increase the block cell counts while preserving dimensions
and grading:

```bash
for refinement in 1 1.25 1.5; do
    for chamber in a b; do
        case="chamber-$chamber-refinement-$refinement"
        python3 examples/run.py "chamber-$chamber" "runs/$case" \
            --refinement "$refinement" --salt-imbalance 0.001 --ranks 8 --results "runs/results/$case"
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
These details are not uniquely specified by the published drawings. The
[original toolkit generator](https://github.com/MathiasGruber/MembraneSimKit/blob/2c4344144dadc3f2b486432724468a91652bbe3a/src/meshLibrary/setupScripts/blockCF042/meshSetup.py)
provides a separate construction with configurable corner radius, distributor
width and inlet height. Its default inlet diameter is 4.5 mm, whereas the paper
specifies 5.5 mm; toolkit defaults alone therefore do not identify the paper's
case settings. This example remains a reconstruction, and exact identity with
the original mesh is unverified.

Default tangential spacing is 0.8 mm with 24 cells through each channel and
wall-normal expansion ratio 100, giving approximately 4 µm at the membrane.
Use `--spacing 0.0004 --layers 48` for a finer grid. Every mesh must pass
`checkMesh`. Solver stopping requires salt residual below 1e-9 and salt
imbalance below 0.1% of membrane transfer. The plotter includes new CFD points
only after convergence and conservation checks; original CFD points are always
labelled separately.

Converged reconstructed CF042 results:

| Flow (mL/min) | Cells | Water, kg/(m² h) | Salt, g/(m² h) | Difference from archived water flux | Salt imbalance / transfer |
|---:|---:|---:|---:|---:|---:|
| 5 | 158,912 | 8.382253 | 5.392870 | −4.7% | 0.0114% |
| 50 | 158,912 | 11.180289 | 7.196734 | −1.7% | 0.0046% |
| 50 | 396,340 | 11.184071 | 7.199158 | −1.6% | <0.1000% |

Relative mass imbalance is below 4e-12. The two 50 mL/min grids differ by 0.034%
in water flux and 0.081% in the draw-surface concentration integral. A third
grid is needed to establish chamber grid independence; the difference from the
archived CFD results remains unresolved.

The original toolkit construction, with the paper's 5.5 mm port and source
wall-normal expansion ratio 5, also converges at 50 mL/min: 826,280 cells,
water 11.188519 kg/(m² h), salt 7.202058 g/(m² h), relative salt imbalance
0.09983% and relative mass imbalance below 1e-12. Its draw-side membrane area
is 0.001653579 m² and its concentration integral is 6.813829e-5 m². Water flux
remains 1.6% below the archived 11.37. This independently constructed mesh is
within 0.040% of the finest accepted reconstruction, but the meshes have
different grading and topology; that agreement is not a grid-convergence
study. The source grading gives a nominal first-cell thickness of 35.90 µm.
This first-cell thickness is a mesh input, not proof of convergence.
Original-construction markers are labelled separately in the figure.

With the same original geometry and expansion ratio 100, the 826,280-cell
baseline converges at water 11.160097 kg/(m² h) and salt 7.183743 g/(m² h), with
relative salt imbalance 0.09911% and relative mass imbalance below 3e-13. Its
nominal first membrane cell is 3.9611 µm. Changing grading from 5 to 100 changes
water flux by −0.254% and the draw-surface concentration integral by −0.347%.
The grading-100 baseline remains 1.85% below the archived water flux.
It is the coarse member of the consistent refinement sequence below.

At expansion ratio 100 and refinement 1.25, the 1,541,448-cell mesh converges
in 6,758 iterations: water 11.166321 kg/(m² h), salt 7.187755 g/(m² h),
relative salt imbalance 0.09986% and relative mass imbalance below 7e-13.
The nominal first membrane cell is 3.2295 µm. Compared with the same-grading
baseline, water changes +0.0558% and the draw-side concentration integral
changes +0.0806% to 6.795666e-5 m². Both paired membrane sides have area
0.001653597 m². This result is 1.79% below the archived water flux.

At refinement 2, the 6,610,240-cell mesh converges in 8,131 iterations:
water 11.173889 kg/(m² h), salt 7.192635 g/(m² h), relative salt imbalance
0.09168% and relative mass imbalance below 9e-15. Its nominal first membrane
cell is 2.0368 µm. From refinement 1.25, water changes +0.06777% and the
draw-side concentration integral changes +0.09559% to 6.8021625e-5 m².
Both membrane sides have area 0.0016536226 m²; independent integration of
60,288 exported faces per side confirms these values.

Using effective spacing proportional to the inverse cube root of cell count,
the three nonuniform grids give apparent orders 1.972 for water and 2.048
for the draw integral. Richardson extrapolation gives water 11.178609
kg/(m² h); the conventional fine-grid GCI estimates, with safety factor 1.25,
are 0.0528% for water and 0.0701% for the draw integral. These are estimates
conditional on an asymptotic single-power error model, not independently
validated uncertainty bounds. Block-count rounding and grading require that
qualification. The finest computed water flux remains 1.725% below the
archived 11.37, substantially larger than the estimated discretization error;
this sequence does not close the historical-input or model discrepancy.

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

### Original toolkit construction

Use `--cf042-geometry original` to run the recovered `blockCF042` construction:

```bash
python3 examples/run.py cf042 runs/cf042-original-50 --cf042-geometry original \
    --flow 50 --salt-imbalance 0.001 --ranks 8 \
    --results runs/results/cf042-original-50
python3 examples/run.py cf042 runs/cf042-original-50-fine --cf042-geometry original \
    --refinement 1.5 --initialize-from runs/cf042-original-50 --ranks 8 \
    --salt-imbalance 0.001 --results runs/results/cf042-original-50-fine
```

To repeat the expansion-100 refinement sequence with the same steady gates:

```bash
python3 examples/run.py cf042 runs/cf042-g100-r1 --cf042-geometry original \
    --wall-normal-expansion 100 --refinement 1 --end 60000 --ranks 4 \
    --salt-imbalance 0.001 --results runs/results/cf042-original-g100-r1
python3 examples/run.py cf042 runs/cf042-g100-r125 --cf042-geometry original \
    --wall-normal-expansion 100 --refinement 1.25 --end 60000 --ranks 4 \
    --initialize-from runs/cf042-g100-r1 --salt-imbalance 0.001 \
    --results runs/results/cf042-original-g100-r125
python3 examples/run.py cf042 runs/cf042-g100-r2 --cf042-geometry original \
    --wall-normal-expansion 100 --refinement 2 --end 60000 --ranks 4 \
    --initialize-from runs/cf042-g100-r125 --salt-imbalance 0.001 \
    --results runs/results/cf042-original-g100-r2
```

The generator preserves the source block topology, a 3 mm corner radius,
29.5 mm distributor width, 0.5 mm inlet buffer and 1 µm merge gap. Override the
first two with `--corner-radius` and `--distributor-width`. Total compartment
height remains `--chamber-height`; the intermediate inlet height is total height
minus channel height. The port diameter defaults to the paper's 5.5 mm, rather
than the toolkit's 4.5 mm. These explicit inputs must still be checked against
the case being reproduced; recovering a generator does not identify every
historical configuration.

The source quarter mesh is merged before mirroring along the length and across
the membrane. The runner splits feed/draw ports, creates paired membrane faces,
and requires `checkMesh` to pass before initializing fields. Identical duplicate
arc declarations are removed for v2606. The default full mesh contains 826,280
cells and two fluid regions. `--refinement` scales the original block counts;
`--spacing` and `--layers` apply only to the other construction. The original
channel grading ratio is 5, and its baseline contains 25 cells through each
channel. Use `--wall-normal-expansion` to test boundary-layer resolution while
preserving the geometry and cell counts. The result records the nominal first
cell height (full first-cell thickness), computed from the geometric progression.
At baseline resolution,
the nominal first cell is 35.90 µm for the source ratio 5 and 3.96 µm for ratio
100. Refinement 2 with ratio 100 gives 50 cells and a nominal 2.04 µm first cell.
These are mesh inputs, not proof of convergence. The source-default mesh needs
its own mesh-convergence study before paper-level conclusions.

Recovering the construction also does not recover every numerical setting.
The toolkit's [steady scheme template](https://github.com/MathiasGruber/MembraneSimKit/blob/2c4344144dadc3f2b486432724468a91652bbe3a/src/templatesFiles/fvSchemes_steadystate)
uses bounded upwind convection for velocity and salt, with cubic interpolation
in momentum diffusion. Its `laplacian(rhoD_AB,m_A)` cubic entry does not match
the coefficient field name `rho*D_AB` in the [archived transport equation](https://github.com/MathiasGruber/MembraneFoam/blob/e16269dccb5b154b2fe22ba61f726c35721050fb/src/solvers/simpleSaltTransport/m_AEqn.H);
that equation therefore selects the template's default linear diffusion.
The maintained examples instead use linear-upwind velocity convection and
van Leer salt convection. Numerical-scheme sensitivity must be distinguished
from geometry and material-parameter sensitivity; these archived templates
alone do not establish the settings of each published case.

Both constructions use the same membrane equations, material properties and
whole-module flow convention. Warm starts require matching construction and
geometry parameters. Collected results record the source revision and inputs;
the dimension plot distinguishes original-toolkit results with separate open
markers and the grading ratio in the legend. Reconstructed geometry uses open
diamonds. Grid selection remains separate for each construction and grading ratio.
Runs with nondefault corner radius or distributor width are excluded from that
reference plot. Mesh quality and startup checks do not establish a converged paper flux or
grid independence.

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
The 13-inlet solution reached its 120,000-step limit without meeting the SIMPLE residual gate: the salt-equation residual stagnated at about 5e-9, five times the 1e-9 tolerance, oscillating without decay over the last 20,000 outer iterations while the velocity and pressure residuals remained below their gates. At the limit the field is steady (water flux 11.747419 kg/(m² h), 1.5% below the archived 11.92 and on the trend between the 11- and 15-inlet points; relative salt imbalance 0.0164% of the membrane transfer; relative mass imbalance 6.9e-13). It is therefore excluded from the accepted table, and the archived 13-inlet point is not reproduced.
At 5 mL/min, the accepted three-inlet, 45° solution gives a water flux of
9.153626 kg/(m² h) on the 1,066,000-cell grid (salt imbalance 0.09956% of
membrane transfer, relative mass imbalance 3.0e-13). The archived 5 mL/min
reference values at 40° and 50° are 9.47 and 9.46 kg/(m² h), so this point
lies 3.3% below the archived reference level, a larger offset than the
50 mL/min points above.
The 10° solution reached the time limit without meeting the 0.1% salt-conservation gate: the relative salt imbalance decayed to 0.115% of the membrane transfer by about a fifth of the run and then plateaued for the remainder, ending at 0.1145%. At the limit the water flux is 12.197 kg/(m² h), 2.4% above the archived 11.91, and the relative mass imbalance is below 3e-13. It is therefore excluded from the accepted table, and the archived 10° point is not reproduced. The 0° solution is numerically unstable: every steady-state trial terminates during the early transient with nonphysical values (negative salt mass fraction on the membrane), and low-order discretization variants remain far above the 0.1% salt-conservation gate when their iteration budget ends. The archived 0° point is therefore not reproduced.
These points do not establish grid independence or reproduce
the complete inlet studies. At 50 mL/min and 45°, the plotter also writes
membrane-field maps for each accepted inlet count.
The 3dBlock runner uses a linear pressure tolerance of 1e-11 and relative
tolerance of 0.001. This avoids driving GAMG below the attainable residual on
the graded inlet meshes. Residual convergence and the 0.1% salt-conservation
requirement must both pass before a result is accepted.

### Spacers

`--spacers` adds an even number of transverse spacers, uniformly spaced along
each channel. Their mesh segments lie at longitudinal positions `length * i / (count + 1)`.
The default `--spacer-placement middle` centres the shape between the membrane
and opposite wall. `up` places it against the outer wall, `down` against the
membrane, and `changing` alternates Down/Up within the length-quarter, then
mirrors that sequence. For six spacers the full order is Down/Up/Down/Down/Up/Down;
it is not strict full-length alternation. At two spacers, Changing and Down
produce the same arrangement. The spacers span
the full width. `--spacer-radius` is an explicit reconstruction parameter. The
[original toolkit configuration](https://github.com/MathiasGruber/MembraneSimKit/blob/2c4344144dadc3f2b486432724468a91652bbe3a/src/meshLibrary/setupScripts/3dBlock/configFile.py)
defines spacer cross-sectional area as `spacerVolume * height²`, with default
`spacerVolume=0.2`. At a 2 mm channel height, this gives a middle-cylinder radius
of 0.504627 mm. `--spacer-area-fraction 0.2` selects that area directly;
`--spacer-shape triangle` uses an equilateral triangle with the same area.
Its vertical base faces increasing x and its centroid lies at mid-height,
matching the source middle-triangle orientation. The defaults provide a
reproducible alternative to the 0.5 mm example below, but are not evidence
of the paper-specific configuration.

```bash
python3 examples/run.py 3dblock runs/block-2-spacers --spacers 2 \
    --spacer-radius 0.0005 --ranks 8
```

For this example radius of 0.5 mm, each cylinder displaces
`π * radius² * width` of fluid per compartment. The mesh grades toward both
the cylinder and channel walls. Use `--refinement` to test spatial convergence
and vary the radius to assess geometric sensitivity. These are reconstructed
middle-spacer cases. The equilateral triangle retains its straight sides and
sharp corners in a conformal mesh; its side length is
`sqrt(4 * area_fraction * height² / sqrt(3))`. This changes the spacer geometry
without importing the legacy inlet mesh or its cell grading. Cylinder and
triangle comparisons should use the same explicit cross-sectional area and
independently converged grids. Wall-contact geometry follows the recovered source boundary construction,
with conformal transitions matching the neighbouring channel cells. The
truncated wall cylinder retains area `(3π/4 + 1/2) * radius²`; an equal-area
wall cylinder therefore has a different radius from a middle cylinder.
Agreement with the published spacer curves remains unverified.

The [2016 paper, §2.3](https://doi.org/10.1016/j.seppur.2015.12.017), describes
approximately 6 million and 26 million cells for its three-inlet, 18-cylinder
3dBlock spacer test. Its 2 µm distance from membrane faces to nearest grid
points refers to that fine spacer mesh, not a specified CF042 mesh resolution.
The nominal first-cell thickness recorded by the CF042 generator is a different
metric and should not be equated with that distance. Reproducing the reported
fine-mesh test requires the corresponding spacer geometry, measured resolution
and transient comparison; cell count alone does not establish spatial or
temporal convergence.

```bash
python3 examples/run.py 3dblock runs/block-cylinder-area --spacers 2 \
    --spacer-area-fraction 0.2 --setup-only
python3 examples/run.py 3dblock runs/block-triangle-area --spacers 2 \
    --spacer-shape triangle --spacer-area-fraction 0.2 --setup-only
```

The standard mesh check passes for the two-triangle setup at area fraction 0.2.
Additional `checkMesh -allTopology -allGeometry` checks flag cell determinants,
concavity and interpolation weights in both the cylinder reconstruction and
triangle mesh. The strongly graded channel and inlet cells therefore require
further mesh-quality and grid-sensitivity assessment; passing the standard
check alone does not establish suitability for a converged paper comparison.
A 20-iteration triangle startup completes with finite, physical salt mass
fractions, but remains far from conservation and residual convergence. No
spacer flux is included in the accepted results.

For wall-contact generation, specify placement explicitly:

```bash
python3 examples/run.py 3dblock runs/block-wall-triangles --spacers 6 \
    --spacer-shape triangle --spacer-placement changing \
    --spacer-area-fraction 0.2 --setup-only
```

All six cylinder/triangle × Up/Down/Changing variants at six spacers pass
the standard mesh check, flow initialization and a bounded ten-iteration
startup. Independently measured areas of both membrane sides match the
expected active area. Salt mass fractions remain physical, but the startups
are far from conservation and residual convergence. Extra geometry diagnostics
retain determinant warnings and the inlet baseline's 660 concave cells; these
checks do not qualify a converged spacer flux or establish grid independence.

Membrane-contact spacers remove permeating membrane area. The historical
boundary conditions normalize transfer by the active paired membrane patch
area, and the toolkit extracts that reported flux directly. At area fraction
0.2, a wall triangle contacts 1.359235 mm along the membrane; a corrected wall
cylinder contacts 0.748456 mm. Six Down triangles reduce nominal active area
by 10.19%; eighteen reduce it by 30.58%. Case metadata records the full placement
order, contact count, nominal area and expected active area. Check these against
actual surface areas before comparison. Active-area flux and total module
transfer are different quantities; a larger reported flux need not mean a
larger total transfer. This source normalization does not identify the inputs
or normalization of every published curve.

The [archived FO water-flux diagnostic](https://github.com/MathiasGruber/MembraneFoam/blob/e16269dccb5b154b2fe22ba61f726c35721050fb/src/boundaryConditions/FO_BC/explicitFOmembraneVelocity/explicitFOmembraneVelocityFvPatchVectorField.C)
uses local `sum` reductions. Under MPI, each reported value is a processor-local
patch average; a global mean requires area-weighted reduction across all ranks.
The toolkit's logfile extractor reads the reported scalar without that reduction.
As an illustration, reconstructing mass flux from the accepted original CF042
expansion-100, refinement-1.25 surface fields gives local means of
9.888–12.192 kg/(m² h), while the global area-weighted mean is 11.166.
This illustrates the reporting issue on current fields, not the historical
execution. Maintained diagnostics reduce globally and use the discrete mass
flux `phi`. The archived reduction is a source defect; its effect on particular
published parallel results remains unverified and does not explain a serial
CF042 discrepancy.

The [source wall-cylinder implementation](https://github.com/MathiasGruber/MembraneSimKit/blob/2c4344144dadc3f2b486432724468a91652bbe3a/src/meshLibrary/setupScripts/3dBlock/spacers/Sphere.py)
uses `0.75 * 3.14159 + 1/2` to resize a cylinder truncated at 45°. In Python 2,
`1/2` evaluates to zero; the resulting retained area exceeds the stated
cross-sectional target by about 21.2%. This identifies a discrepancy in the
archived source. It does not establish which revision or inputs generated the
paper's wall-spacer fluxes, or prove those fluxes erroneous.

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
