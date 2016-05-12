#!/bin/sh

# Load OpenFOAM
# Older OpenFOAM version is used for setting up mesh here, but newer should work also
module load OpenFoam/2.3.1/gcc-4.8.4-openmpi

# Remove previous 0 dir
rm -rf 0

# First run blockMesh to get the basic mesh
blockMesh

# Move the mesh to the right position
transformPoints -translate "(-0.04 0 0)"

# Load/run the first mirrorMeshDict
cp system/mirrorMeshDict_X system/mirrorMeshDict
mirrorMesh

# Load/run the second mirrorMeshDict
cp system/mirrorMeshDict_Y system/mirrorMeshDict
mirrorMesh

# Load a fixed boundary file. On each mirror mesh, the inlet/outlet patches are doubled
# in size, first in the X-direction and then Y-direction. These are then split manually (or via script)
# to have the right names (outletTop,outletBottom, inletTop, inletBottom)
cp constant/polyMesh/boundaryFixed constant/polyMesh/boundary

# Create a baffle on the membrane
topoSet && setsToZones -noFlipMap && createBaffles

# Move time directory from 1 to 0
mv 1/ 0/

# Copy in the initial fields
cp initialFields/* 0/

# Run the setFields dictionary for better initial fields
setFields