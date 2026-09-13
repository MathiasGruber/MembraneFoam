# Native GPU execution

Native offloading requires the official [OpenFOAM ECSE](https://gitlab.com/openfoam/gpu/openfoam-ecse)
`feature-gpu` source and NVIDIA HPC SDK. A release CPU binary does not offload.
The API-2606 revision used for compatibility checks is
`f36870083be0f391cf44c4aa6e1623214282e271`.
Follow the upstream [build guide](https://gitlab.com/openfoam/gpu/openfoam-ecse/-/wikis/build)
and [dependency guide](https://gitlab.com/openfoam/gpu/openfoam-ecse/-/wikis/dependencies)
for the target architecture, compiler, and memory model. OpenFOAM GPU toolchain
installation and architecture-specific patches are outside this library.

Build MembraneFoam in the GPU environment, keeping CPU and GPU libraries separate:

```bash
source /path/to/OpenFOAM-gpu/etc/bashrc WM_COMPILER=Nvidia-gpu
./Allwmake
```

Prepare a case with CPU mesh utilities using `examples/run.py --setup-only`,
then run `simpleSaltTransport -case /path/to/case` in the GPU environment.
For transient cases, include `--transient` and `--time-step` during preparation,
then run `pisoSaltTransport -case /path/to/case` in the GPU environment.
Use `examples/run.py --collect-only` to validate and export the solved case.
The supported configuration uses one GPU per simulation. Membrane pairing,
root solving, and some linear-algebra stages remain on the CPU.

Compare identical initial fields, dictionaries, convergence criteria, and integral
fluxes when evaluating CPU/GPU performance. Exclude compilation and meshing from
solver timings and use warm-up runs and repeated measurements. GPU execution
does not imply a speedup, especially for small meshes.
