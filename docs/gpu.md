# Running on a GPU

Use native OpenFOAM ECSE offload with AmgX for pressure, `dFFI` inner precision,
hierarchy reuse and exact-matrix reuse. Keep momentum and salt on native
PBiCGStab with diagonal preconditioning. Initialize with 500 consistent SIMPLE
iterations, then restore the intended equations and converge them. This is the
best validated converged chamber configuration described below.

OpenFOAM fields, assembly, membrane physics and final residual acceptance remain
double precision. Float arithmetic is confined to AmgX pressure corrections.
The adapter supports one GPU in one serial process, a static mesh and 32-bit
labels. Pressure must be positive definite after sign normalization. Membrane
pairing and nonlinear boundary calculations include CPU work.

## Build

Use [OpenFOAM ECSE](https://gitlab.com/openfoam/gpu/openfoam-ecse), branch
`feature-gpu`, with NVIDIA HPC SDK. The compatible API-2606 revision is
`f36870083be0f391cf44c4aa6e1623214282e271`. Follow the upstream
[build guide](https://gitlab.com/openfoam/gpu/openfoam-ecse/-/wikis/build) and
[dependency guide](https://gitlab.com/openfoam/gpu/openfoam-ecse/-/wikis/dependencies).

The GB10 configuration uses NVIDIA HPC SDK 26.5, double fields and 32-bit labels.
Set these definitions in the separate OpenFOAM checkout's
`wmake/rules/General/Nvidia-gpu/c++` and `link-c++` files:

```make
CC := nvc++$(COMPILER_VERSION) -std=c++17 -stdpar=gpu -gpu=cc121,mem:unified:nomanagedalloc -DFOAM_OFFLOAD
LINKLIBSO = $(CC) $(c++FLAGS) -shared -fPIC -nonvmalloc
LINKEXE = $(CC) $(c++FLAGS) -nonvmalloc -Xlinker --add-needed
```

At this revision, include `fvMesh.H` before the scheme header in
`src/finiteVolume/finiteVolume/snGradSchemes/orthogonalSnGrad/orthogonalSnGrads.C`
and `src/finiteVolume/interpolation/surfaceInterpolation/schemes/localMin/localMin.C`.
Set `UMPIRE_ARCH_PATH` to the installed Umpire dependency and use
`FOAM_MEMORY_POOL="managed; size=1024; incr=1024"`. The explicit memory pool and
compiler system-allocation mode govern different allocations. See the
[NVIDIA memory-mode reference](https://docs.nvidia.com/hpc-sdk/archive/26.5/compilers/hpc-compilers-user-guide/index.html).

Build AmgX using the [pinned dependency instructions](amgx.md#build). The GB10
dependency uses CUDA 13.0.88, GCC 13.3 and architecture 121; the OpenFOAM adapter
uses the HPC SDK's CUDA 13.2 headers/runtime. From the MembraneFoam repository
root, build the application and optional plugin:

```bash
source /path/to/OpenFOAM-gpu/etc/bashrc WM_COMPILER=Nvidia-gpu
export FOAM_MEMORY_POOL="managed; size=1024; incr=1024"
export AMGX_SOURCE=/path/to/AMGX
export AMGX_BUILD=/path/to/directory/containing/libamgxsh.so
export CUDA_HOME=/path/to/compiler-compatible/cuda
./Allwmake
./src/linearSolvers/amgx/Allwmake
```

Keep CPU and GPU installations and output libraries separate. Use CPU mesh
utilities and the CPU membrane library for preparation. Verify loaded libraries
and run the [operator checks](amgx.md#operator-verification) in the GPU environment
before using a different compiler/CUDA combination.

## Prepare and configure the case

In the CPU OpenFOAM environment, run these commands from the MembraneFoam
repository root with a new absolute case path:

```bash
case_dir=/absolute/path/to/chamber
python3 examples/run.py chamber-a "$case_dir" --setup-only \
    --end 20600 --salt-imbalance 0.001
cp src/linearSolvers/amgx/configs/pcg-aggregation-reuse.json \
    "$case_dir/system/amgx.json"
```

Set the GPU solvers and chamber accuracy requirements explicitly:

```bash
foamDictionary "$case_dir/system/controlDict" -entry libs \
    -set '("libDHIBoundaryConditions.so" "libmembraneAmgX.so")'
foamDictionary "$case_dir/system/controlDict" -entry writeFormat -set binary
foamDictionary "$case_dir/system/fvSolution" -entry solvers -set '
{
 p
 {
  solver AmgX; tolerance 1e-12; relTol 0.001; maxIter 2000;
  AmgX
  {
   configFile "$FOAM_CASE/system/amgx.json";
   mode dFFI; verifyCsr true;
   reuseHierarchy true; reuseUnchangedMatrix true;
  }
 }
 pFinal { $p; tolerance 1e-10; relTol 0; }
 "(U|UFinal|m_A)"
 {
  solver PBiCGStab; preconditioner diagonal;
  tolerance 1e-15; relTol 0; maxIter 2000;
 }
}'
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.consistent -set false
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.consistentAudit -set false
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.residualControl \
    -set '{ p 1e-6; U 1e-8; m_A 1e-9; }'
foamDictionary "$case_dir/system/fvSolution" \
    -entry SIMPLE.maxRelativeSaltImbalance -set 0.001
foamDictionary "$case_dir/system/fvSolution" \
    -entry SIMPLE.validationConsecutiveChecks -set 100
cp "$case_dir/system/fvSolution" "$case_dir/system/fvSolution.target"
cp "$case_dir/system/fvSchemes" "$case_dir/system/fvSchemes.target"
```

For another case, preserve its required pressure, momentum, salt and outer
accuracy. PCG needs a suitable symmetric positive-definite pressure operator.
Use the [consistent method's supported discretizations](consistent-initialization.md#supported-discretization).
The chamber has two non-orthogonal correctors and `Gauss linear corrected`
pressure Laplacians.

## Run the steady solve

Switch to the GPU environment used to build the application. Keep `case_dir`
set to the prepared case. Use one process and one OpenMP/BLAS thread:

```bash
source /path/to/OpenFOAM-gpu/etc/bashrc WM_COMPILER=Nvidia-gpu
export FOAM_MEMORY_POOL="managed; size=1024; incr=1024"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
```

Initialize for 500 iterations with consistent pressure coupling, upwind salt
convection and relaxation factors `p=1`, `U=0.7`, `m_A=0.6`:

```bash
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.consistent -set true
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.residualControl -set '{}'
foamDictionary "$case_dir/system/fvSolution" -entry relaxationFactors \
    -set '{ fields { p 1; } equations { U 0.7; m_A 0.6; } }'
foamDictionary "$case_dir/system/fvSchemes" \
    -entry 'divSchemes.div(phi,m_A)' -set 'Gauss upwind'
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 500
foamDictionary "$case_dir/system/controlDict" -entry writeInterval -set 500
simpleSaltTransport -case "$case_dir" > "$case_dir/log.startup" 2>&1
```

Restore the intended equations, relaxation and linear targets. Run 100 iterations
of these equations before enabling final convergence checks:

```bash
cp "$case_dir/system/fvSolution.target" "$case_dir/system/fvSolution"
cp "$case_dir/system/fvSchemes.target" "$case_dir/system/fvSchemes"
foamDictionary "$case_dir/system/fvSolution" -entry SIMPLE.residualControl -set '{}'
foamDictionary "$case_dir/system/controlDict" -entry startFrom -set latestTime
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 600
foamDictionary "$case_dir/system/controlDict" -entry writeInterval -set 100
simpleSaltTransport -case "$case_dir" > "$case_dir/log.probe" 2>&1

cp "$case_dir/system/fvSolution.target" "$case_dir/system/fvSolution"
foamDictionary "$case_dir/system/controlDict" -entry endTime -set 20600
foamDictionary "$case_dir/system/controlDict" -entry writeInterval -set 1000
simpleSaltTransport -case "$case_dir" > "$case_dir/log.solver" 2>&1
```

Require steady convergence before collecting the result in the CPU environment:

```bash
python3 examples/run.py chamber-a "$case_dir" --collect-only \
    --end 20600 --salt-imbalance 0.001
```

Final acceptance also requires physical boundary mass imbalance at most `1e-8`,
relative water-flux span at most `1e-4` over the final 100 iterations and bulk
salt bounds `[-1e-7, 1 + 1e-7]`. Paired CPU/GPU verification requires global
unweighted field L2 differences below `1e-7`, maximum velocity-vector difference
below `1e-8 m/s` and membrane water means within `0.1%`. Consecutive solver
checks do not replace these independent field, flux and physical-bounds checks.

## Validated performance and limits

For the 326,720-cell chamber with identical generated initial fields and the
method above:

| Backend | Time to accepted solution |
|---|---:|
| CPU, 20 MPI ranks | 79.49 min |
| GB10 GPU, AmgX pressure | 34.65 min |

The paired speedup is `2.29×`. The CPU reference uses GCC 13.3,
GAMG/GaussSeidel pressure with default coarsest-level settings and native
PBiCGStab/DILU for momentum and salt. Its linear targets, iteration caps and
final accuracy match the GPU case. MPI uses a simple `(4 5 1)` decomposition
with `preserveBaffles`, `--bind-to none` and one OpenMP/BLAS thread per rank.

Timing includes case preparation/configuration, CPU decomposition, startup,
handoff, the probe and final convergence. Compilation, meshing and shared initial
potential-flow construction are excluded. Maximum CPU/GPU velocity difference
is below `5.6e-9 m/s`, the largest relative field L2 difference is below `1.18e-8`
and raw physical mass imbalance is below `5.7e-14`.

This result is one paired sample on one mesh. It does not establish a speedup
for every case, full GPU utilization or uniqueness of the steady solution.
Keep identical initial fields and final accuracy when comparing configurations,
include initialization cost and use repeated unprofiled converged timings to
estimate variation.
