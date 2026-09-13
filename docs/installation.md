# Installation on Linux

Use **OpenCFD OpenFOAM v2606**, double precision and 32-bit mesh labels.
Source packages are available from the [official release downloads](https://dl.openfoam.com/source/v2606/).

For Ubuntu 24.04, install the build tools and libraries:

```bash
sudo apt-get update
sudo apt-get install build-essential flex bison libfl-dev libopenmpi-dev \
    zlib1g-dev libreadline-dev libncurses-dev libxt-dev libscotch-dev
```

Download and unpack the archive into your chosen build directory, then open a
Bash shell and source its environment. Source OpenFOAM **before** enabling
shell `set -e`; its environment scripts are not designed for that mode.

```bash
source /path/to/OpenFOAM-v2606/etc/bashrc
cd /path/to/MembraneFoam
./Allwmake
tests/run_openfoam.sh
```

Install OpenFOAM using its distribution packages or upstream source-build
instructions before building MembraneFoam. Keep libraries from different
OpenFOAM distributions and versions in separate environments.

The chamber examples require `blockMesh`, `mirrorMesh`, `topoSet`,
`createBaffles`, `createPatch`, and `checkMesh`. Refined 2012 chambers also use
`collapseEdges`; include this utility when building a reduced OpenFOAM toolset.

Every shell that runs a solver needs the selected OpenFOAM environment. Executables
and the membrane library install into `FOAM_USER_APPBIN` and `FOAM_USER_LIBBIN`.
Load `libDHIBoundaryConditions.so` in the case's `controlDict`; generated examples
already do this. Use a fresh generated case rather than copying a legacy mesh
and assuming its dimensions, patch names and boundary-condition settings match.

For native GPU offloading, follow [the separate GPU build](gpu.md). A CPU binary
continues to use the CPU when run on a machine with a GPU.
