#!/bin/bash
set -euo pipefail

example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

repo_root=${SCALEUPROM_ROOT:-}
if [[ -z ${repo_root} ]]; then
    search_dir=${example_dir}
    while [[ ${search_dir} != / ]]; do
        if [[ -d ${search_dir}/utils/python && -d ${search_dir}/examples/stokes/meshes ]]; then
            repo_root=${search_dir}
            break
        fi
        search_dir=$(dirname -- "${search_dir}")
    done
fi

if [[ -z ${repo_root} ]]; then
    echo "Could not locate the scaleupROM source root." >&2
    exit 1
fi

converter=${SCALEUPROM_GMSH2MFEM:-${repo_root}/build/utils/gmsh2mfem}
if [[ ! -x ${converter} ]]; then
    echo "Missing ${converter}; build scaleupROM before preparing this example." >&2
    exit 1
fi

cd "${example_dir}"
mkdir -p meshes configs

for mesh in square-circle.msh.mfem square-square.msh.mfem \
            square-triangle.msh.mfem square-star.msh.mfem; do
    cp "${repo_root}/examples/stokes/meshes/${mesh}" "meshes/${mesh}"
done

# Promote the linear triangular empty-cell mesh to the exact cubic mesh used
# by the paper. gmsh2mfem appends ".mfem" to its input name.
cp "${repo_root}/test/meshes/square.tri.mesh" meshes/square.tri.mesh
"${converter}" -m meshes/square.tri.mesh -o 3
mv meshes/square.tri.mesh.mfem meshes/square.tri-elem.mesh

python3 generate_configs.py --output configs
