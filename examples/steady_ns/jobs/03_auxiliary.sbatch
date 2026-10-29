#!/bin/bash
# Manual recovery only: train_rom normally calls aux_train_rom itself.
#SBATCH --job-name=ns-supremizer
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
    "main/mode=aux_train_rom:main/use_rom=true:basis/number_of_basis=40:basis/number_of_supremizer=40"
