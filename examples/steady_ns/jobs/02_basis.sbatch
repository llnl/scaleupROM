#!/bin/bash
#SBATCH --job-name=ns-basis
#SBATCH --nodes=5
#SBATCH --ntasks=36
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --partition=pbatch
#SBATCH --account=genddfem
#SBATCH --output=logs/%x-%j.out

# shellcheck source=jobs/common.bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.bash"

srun --ntasks=36 "${MAIN_BIN}" -i steady_ns.yml -f \
    "main/mode=train_rom:main/use_rom=true"
