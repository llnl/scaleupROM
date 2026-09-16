#!/bin/bash
# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT
#
# Generates the input meshes for the high-Re channel example and creates the
# output directories the solver writes into.
#
# Note: assumes this is run from the build/examples/highReChannel/ directory.
#
#   ./setup_highrechannel.sh                 # both default resolutions
#   ./setup_highrechannel.sh 0.015625        # an extra, finer mesh
#
# Must be run once before any of the bump2d.*.yml configs will load: the configs
# reference .msh.mfem files that this script produces.

set -euo pipefail

GMSH2MFEM=../../utils/gmsh2mfem

if [[ ! -x ${GMSH2MFEM} ]]; then
    echo "error: ${GMSH2MFEM} not found."
    echo "CMake only builds utils/gmsh2mfem when find_program(GMSH gmsh) succeeds,"
    echo "so make sure gmsh is on PATH when configuring, then rebuild."
    exit 1
fi

if ! python3 -c "import gmsh" 2>/dev/null; then
    echo "error: the gmsh python module is not importable."
    echo "Install it with 'pip3 install gmsh' (the apt 'gmsh' package ships the"
    echo "CLI only, which is enough for CMake but not for generate_mesh.py)."
    exit 1
fi

# h = 0.0625 -> Re1e5 and the smoke test; h = 0.03125 -> Re500.
# See README.md for why the lower-Re case gets the finer mesh.
resolutions=("$@")
if [[ ${#resolutions[@]} -eq 0 ]]; then
    resolutions=(0.0625 0.03125)
fi

for h in "${resolutions[@]}"; do
    msh="meshes/bump2d.h${h}.msh"
    python3 generate_mesh.py --h "${h}" --out "${msh}" --check

    # -o 1 keeps the mesh straight-sided. gmsh2mfem defaults to order 3, which
    # would promote the nodes to a discontinuous cubic space for no benefit on a
    # rectilinear geometry.
    ${GMSH2MFEM} -m "${msh}" -o 1 > /dev/null
    echo "  converted      : ${msh}.mfem"
done

# SaveSolutionWithTime calls H5Fcreate without creating the directory first.
mkdir -p paraview restart logs

echo
echo "Setup complete. Next: ../../bin/main -i bump2d.smoke.yml"
