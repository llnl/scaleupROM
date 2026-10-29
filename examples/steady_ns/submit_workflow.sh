#!/bin/bash
set -euo pipefail

example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "${example_dir}"
mkdir -p logs

extra_sbatch_args=("$@")

submit() {
    local dependency=$1
    local script=$2
    local result
    local args=(--parsable)
    if [[ -n ${dependency} ]]; then
        args+=("--dependency=afterok:${dependency}")
    fi
    result=$(sbatch "${args[@]}" "${extra_sbatch_args[@]}" "${script}")
    printf '%s' "${result%%;*}"
}

mesh_id=$(submit "" jobs/00_mesh.sbatch)
samples_id=$(submit "${mesh_id}" jobs/01_samples.sbatch)
basis_id=$(submit "${samples_id}" jobs/02_basis.sbatch)
build_id=$(submit "${basis_id}" jobs/04_build_rom.sbatch)
fom_id=$(submit "${mesh_id}" jobs/05_fom_16x16.sbatch)
rom_id=$(submit "${build_id}:${fom_id}" jobs/06_rom_compare_16x16.sbatch)

printf 'mesh:       %s\n' "${mesh_id}"
printf 'samples:    %s\n' "${samples_id}"
printf 'basis:      %s\n' "${basis_id}"
printf 'ROM build:  %s\n' "${build_id}"
printf '16x16 FOM:  %s\n' "${fom_id}"
printf '16x16 ROM:  %s\n' "${rom_id}"
