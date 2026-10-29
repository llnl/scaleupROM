#!/bin/bash
#SBATCH --job-name=ns-mesh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --partition=pbatch
#SBATCH --account=genddfem
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
example_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "${example_dir}"
mkdir -p logs
./prepare_meshes.sh

