#!/bin/bash
#SBATCH --job-name=ns-fom-16
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

mkdir -p results/paraview/fom
"${MAIN_BIN}" -i steady_ns.yml -o results/fom_timings.h5 -f \
    "main/mode=single_run:main/use_rom=false:mesh/component-wise/global_config=configs/paper.16x16.h5:model_reduction/compare_solution/enabled=false:save_solution/enabled=true:save_solution/file_path/prefix=results/fom_solution:visualization/enabled=true:visualization/file_path/prefix=results/paraview/fom/steady_ns"
