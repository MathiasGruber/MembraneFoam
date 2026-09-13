# 2012 membrane chambers

From a Bash shell with OpenFOAM v2606 and MembraneFoam built:

```bash
python3 examples/run.py chamber-b runs/chamber-b --ranks 8 \
    --results runs/results/chamber-b
```

Run from the repository root. Use `chamber-a` for chamber A. `--setup-only` assembles and
checks the mesh and initializes the flow without launching the transport solver.
Each case directory must be new. The launcher logs generation, mirroring,
baffle creation, port selection, mesh checks and potential-flow initialization.
The CPU MPI decomposition uses `preserveBaffles` to keep paired faces local.

The generators derive from `createSMTCdict.py` (A) and `createFO21Dict.py` (B).
See [geometry and validation details](../../../docs/reproduction.md) for dimensions,
chamber-A geometry adjustments, numerical settings, and comparison results.
