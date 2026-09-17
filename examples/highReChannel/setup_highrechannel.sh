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
#   ./setup_highrechannel.sh                        # everything the configs need
#   ./setup_highrechannel.sh --case corner2d        # just one geometry
#   ./setup_highrechannel.sh -h 0.015625            # an extra, finer resolution
#   ./setup_highrechannel.sh --tags bfs             # skip the channel_flow variants
#
# Must be run once before any of the *.yml configs will load: they reference
# .msh.mfem files that this script produces.
#
# Two tag variants of each mesh are generated because backward_facing_step and
# channel_flow disagree about what boundary tag 2 means -- see README.md.

set -euo pipefail

CASES=(bump2d corner2d)
RESOLUTIONS=(0.0625 0.03125)
TAGS=(bfs chan)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --case)  CASES=("$2"); shift 2 ;;
        --tags)  TAGS=("$2"); shift 2 ;;
        -h|--h)  RESOLUTIONS=("$2"); shift 2 ;;
        --help)  sed -n '6,21p' "$0"; exit 0 ;;
        *)       echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

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

# h = 0.0625 -> the Re=1e5 configs and the smoke tests; h = 0.03125 -> Re=500.
# See README.md for why the lower-Re case gets the finer mesh.
for case in "${CASES[@]}"; do
    for h in "${RESOLUTIONS[@]}"; do
        for tags in "${TAGS[@]}"; do
            msh="meshes/${case}.h${h}.${tags}.msh"
            python3 generate_mesh.py --case "${case}" --tags "${tags}" \
                                     --h "${h}" --out "${msh}" --check \
                | grep -E "^wrote|elements |boundary edges |check |FAIL"

            # -o 1 keeps the mesh straight-sided. gmsh2mfem defaults to order 3,
            # which would promote the nodes to a discontinuous cubic space for no
            # benefit on a rectilinear geometry.
            ${GMSH2MFEM} -m "${msh}" -o 1 > /dev/null
            echo "  converted      : ${msh}.mfem"
        done
    done
done

# SaveSolutionWithTime calls H5Fcreate without creating the directory first.
mkdir -p paraview restart logs

echo
echo "Setup complete. Next: ../../bin/main -i bump2d.smoke.yml"
echo "              or:    ../../bin/main -i corner2d.smoke.yml"
