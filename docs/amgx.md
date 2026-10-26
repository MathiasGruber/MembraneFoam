# AmgX pressure solver

The `libmembraneAmgX` plugin provides an `AmgX` runtime solver for OpenCFD
OpenFOAM v2606. It supports one GPU in one serial process, a static finite-volume
mesh, double-precision OpenFOAM fields and 32-bit labels. Coupled interfaces,
MPI, moving meshes, block matrices and positive `minIter` are unsupported.
The membrane library and default CPU build do not depend on AmgX.

## Build

Build [NVIDIA AmgX](https://github.com/NVIDIA/AMGX) as a shared library for the
target CUDA architecture. The adapter uses the public C API. The reference
AmgX revision is `91a8413ef267b1c32aff4014c02820e1c5897ac2`. Classical AMG with
CUDA 13 requires the kernel-instantiation repair in commit
[`b990b94d38b345ec828601198703af00069db142`](https://github.com/NVIDIA/AMGX/commit/b990b94d38b345ec828601198703af00069db142),
applied to that revision. It adds the integer-sequence device instantiation in
`src/amg.cu`. Pin the base and repair for a reproducible dependency build.
Validate the standalone AmgX examples before building the adapter.

Source the intended OpenFOAM environment, then set dependency locations:

```bash
export AMGX_SOURCE=/path/to/AMGX
export AMGX_BUILD=/path/to/directory/containing/libamgxsh.so
export CUDA_HOME=/path/to/compiler-compatible/cuda
./src/linearSolvers/amgx/Allwmake
```

Keep CPU and GPU OpenFOAM builds and output directories separate. The CUDA
headers and runtime must match the compiler toolchain. A native GPU build uses
NVIDIA stdpar for coefficient conversion and update loops; the allocation
model must make those buffers accessible to the GPU. See the upstream
[GPU build instructions](gpu.md) and NVIDIA's
[stdpar memory documentation](https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/#using-stdpar).
A CPU OpenFOAM build keeps those loops serial and delegates only the sparse
solve to AmgX. Both configurations still require a CUDA-capable GPU at runtime.

## Configure

Add the plugin to the case's existing `libs` list in `system/controlDict`,
retaining the membrane boundary library and any other required libraries:

```foam
libs ("libDHIBoundaryConditions.so" "libmembraneAmgX.so");
```

Copy `src/linearSolvers/amgx/configs/pcg-aggregation-reuse.json` to
`system/amgx.json`. It uses PCG with aggregation AMG and retained hierarchy
topology, requiring a symmetric positive-definite pressure operator after
sign normalization. Use the following pressure configuration while preserving
the case's required `tolerance`, `relTol` and `maxIter`:

```foam
p
{
    solver      AmgX;
    tolerance   1e-12;
    relTol      0.001;
    maxIter     2000;
    AmgX
    {
        configFile "$FOAM_CASE/system/amgx.json";
        mode       dFFI;
        verifyCsr  true;
        reuseHierarchy true;
        reuseUnchangedMatrix true;
    }
}
```

Configure `pFinal` with its separate original target. The complete chamber
configuration and native momentum/salt settings are in the
[GPU execution guide](gpu.md#prepare-and-configure-the-case).
AmgX mode names encode vector precision before matrix precision:

| Mode | Vectors | Matrix coefficients |
|---|---|---|
| `dDDI` | double | double |
| `dDFI` | double | float |
| `dFFI` | float | float |

The OpenFOAM fields, residual calculation and correction accumulation stay in
double precision in all modes. A mode must also be supported by the selected
AmgX algorithms and dependency build. The reference AmgX revision rejects
`dDFI` on CUDA 10.1 and later.
`dFFI` still permits float inner corrections with double OpenFOAM residuals.
Float matrix overflow is an error.

## Method and limitations

The plugin caches the CSR sparsity map, conversion buffers and AmgX handles in
the mesh registry. By default, every solve refreshes coefficients and rebuilds
the numerical AMG hierarchy. Construction uses compact row offsets and one
temporary entry array, with independent values and origin maps for each field.
Repeated columns are summed; those matrices rebuild the conversion on coefficient
updates so that every duplicate contribution is refreshed.

Two optional `AmgX` settings reduce repeated work:

```foam
reuseHierarchy       true;
reuseUnchangedMatrix true;
```

Both default to `false` and can be enabled independently. `reuseHierarchy` uses
AmgX numeric resetup after uploading changed coefficients. The configuration's
`structure_reuse_levels` controls retained AMG topology: `0` rebuilds it, while
`-1` retains all levels. The supplied pressure configuration selects `-1`.
Numerical coarse operators are still refreshed; retaining topology does not
permit applying obsolete matrix coefficients. Changing coefficients can also
make an old hierarchy less effective, so compare iteration counts and elapsed
time with full setup.

`reuseUnchangedMatrix` skips coefficient upload and setup only when every
current double-precision CSR coefficient exactly equals the cached matrix.
It retains an extra double buffer of about eight bytes per nonzero, as well as
comparison/copy work. It never approximates equality with a tolerance or hash.
The RHS is always refreshed, corrections start from zero, and acceptance still
uses the current original OpenFOAM operator. This option is useful only when
enough successive solves share a matrix to offset its cost.

Vector components can have different scalar diagonals because boundary terms
depend on the component. Sharing a matrix across components requires equality
of the assembled operators; identical sparsity alone is insufficient.

Changing the mode, configuration contents or reuse policy replaces the cached
state for that field. Mesh destruction releases its solver, vector and matrix
handles. All cached fields share one AmgX resource owned by the runtime:
resource destruction releases process-wide device pools and library handles,
so it must follow destruction of every dependent field handle. The runtime
remains in the OpenFOAM Time registry across mesh lifetimes and finalizes
before the application unloads libraries.

The runtime retains an independent copy of its first solver configuration for
the shared resource. Keep resource-wide settings, including memory-pool,
stream and determinism options, consistent across configuration files. Solver
algorithms, precision and inner stopping rules remain independent per field.
Changing a field's solver configuration does not recreate the shared resource.
The plugin owns its AmgX runtime;
combining it with another independently managed AmgX integration is unsupported.

AmgX computes an inexact residual correction with a zero initial guess. After
each correction, OpenFOAM recomputes the true residual with its original matrix
and the normalization fixed at the start of that solve. Only the original
OpenFOAM absolute/relative target permits acceptance. Breakdown, stagnation,
the correction limit or the OpenFOAM iteration cap cause an error. AmgX's own
internal norm and tolerance govern inner work and do not replace that target.
With `verifyCsr true`, the first three operators for each cached state are also
compared against OpenFOAM matrix-vector products, including coefficient updates.

This accelerates linear algebra. Finite-volume assembly and custom membrane
work still depend on the chosen OpenFOAM backend. Fixed-iteration throughput,
restart tests and time to a physically accepted solution measure different
things. Compare identical inputs, original convergence and conservation targets,
fields, surface integrals and repeated elapsed times against the fastest tested
CPU configuration; see [performance and accuracy limits](gpu.md#validated-performance-and-limits).

## Operator verification

Build the manufactured-operator regression in the same environment as the plugin:

```bash
wmake tests/amgx
cp -R tests/amgx/case /path/to/new-verification-case
blockMesh -case /path/to/new-verification-case
testMembraneAmgX -case /path/to/new-verification-case \
    -config "$PWD/src/linearSolvers/amgx/configs/fgmres-aggregation.json"
testMembraneAmgX -case /path/to/new-verification-case \
    -config "$PWD/src/linearSolvers/amgx/configs/fgmres-aggregation-reuse.json" \
    -reuseHierarchy -reuseUnchangedMatrix
testMembraneAmgX -case /path/to/new-verification-case \
    -config "$PWD/src/linearSolvers/amgx/configs/fgmres-aggregation-reuse.json" \
    -momentumConfig "$PWD/src/linearSolvers/amgx/configs/pbicgstab-jacobi.json" \
    -simultaneousFields -reuseHierarchy -reuseUnchangedMatrix
```

Prepare the mesh with CPU utilities before switching to a native GPU environment.
The supplied fixture is a 2,048-cell cube with uncoupled walls; another static,
uncoupled mesh with the usual control/discretization dictionaries can be used.
No solution fields are read.
The utility checks known solutions of nonsymmetric, diagonally dominant
operators, alternating changed coefficients and changed RHS, a negative-diagonal
operator, repeated cached updates, double/float vector modes and mesh
teardown/reinitialization. The default modes give 60 known-solution checks.
`-simultaneousFields` gives 720 checks with four cached field states and
overlapping meshes. It destroys the second mesh and solves again on the first
before releasing the shared resource. `-momentumConfig` supplies a separate
solver configuration for the three momentum states. For instrumented runs,
`-compact` reduces this workload to 180 checks while retaining configuration,
precision, sign, matrix/RHS changes and all mesh/state destruction paths.
Test each reuse flag separately as well as together, and retain the full-setup
control. Run the simultaneous-field test under Compute Sanitizer with
`--tool memcheck --leak-check full --error-exitcode 86` and `-compact`; require
a complete test, zero errors and zero leaked bytes. Also check a large-mesh
run across enough steps to refresh coefficients and exercise the selected
setup/resetup path, including normal process exit. Confirm that its diagnostics
show the intended updates: a single step may cover only initial setup. Add
`-mixed` to verify `dDFI` when the AmgX build supports it.
Operator verification does not establish membrane
convergence, physical equivalence or a performance improvement.

## Consecutive steady convergence checks

`simpleSaltTransport` can require consecutive passes of its original outer
residual criteria and integral salt balance:

```foam
SIMPLE
{
    // Retain the case's original residualControl and other SIMPLE settings.
    maxRelativeSaltImbalance 0.001;
    validationConsecutiveChecks 100;
}
```

The default is one check; valid counts are 1 to 100,000. A residual or salt-balance
failure resets the counter. This does not change linear tolerances or replace
independent flux stability, field comparisons and physical-bounds checks.
Small global residuals and stable membrane flux can coexist with localized
velocity differences. Specify the required local-field accuracy separately.

The adapter's derivation and license are recorded in
[`NOTICE`](../src/linearSolvers/amgx/NOTICE) and
[`COPYING`](../src/linearSolvers/amgx/COPYING).
