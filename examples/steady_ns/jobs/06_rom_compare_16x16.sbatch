#!/bin/bash
#SBATCH --job-name=ns-rom-16
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

mkdir -p results/paraview/rom
"${MAIN_BIN}" -i steady_ns.yml -o results/rom_vs_fom.h5 -f \
    "main/mode=single_run:main/use_rom=true:mesh/component-wise/global_config=configs/paper.16x16.h5:basis/number_of_basis=40:basis/number_of_supremizer=40:model_reduction/nonlinear_handling=tensor:model_reduction/save_operator/prefix=rom/steady_ns_nb40:model_reduction/compare_solution/enabled=true:model_reduction/compare_solution/load_solution=true:model_reduction/compare_solution/fom_solution_file=results/fom_solution.h5:save_solution/enabled=false:visualization/enabled=true:visualization/file_path/prefix=results/paraview/rom/steady_ns"

python3 report_results.py
