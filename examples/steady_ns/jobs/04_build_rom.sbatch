#!/bin/bash
#SBATCH --job-name=ns-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --partition=pbatch
#SBATCH --account=genddfem
#SBATCH --output=logs/%x-%j.out

# shellcheck source=jobs/common.bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.bash"

"${MAIN_BIN}" -i steady_ns.yml -f \
    "main/mode=build_rom:main/use_rom=true:basis/number_of_basis=40:basis/number_of_supremizer=40:model_reduction/nonlinear_handling=tensor:model_reduction/save_operator/prefix=rom/steady_ns_nb40"
