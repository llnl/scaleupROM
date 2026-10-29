#!/bin/bash
#SBATCH --job-name=ns-samples
#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --partition=pbatch
#SBATCH --account=genddfem
#SBATCH --output=logs/%x-%A_%a.out

# shellcheck source=jobs/common.bash
source "$(dirname -- "${BASH_SOURCE[0]}")/common.bash"

batch_id=${SLURM_ARRAY_TASK_ID:-${1:-}}
if [[ -z ${batch_id} || ${batch_id} -lt 1 || ${batch_id} -gt 5 ]]; then
    echo "Batch index must be in [1, 5]." >&2
    exit 2
fi

mkdir -p "snapshots/batch-${batch_id}"
# RandomSampleGenerator seeds from wall time. Offset simultaneous array tasks
# so that independent batches cannot receive the same one-second seed.
sleep $((5 * batch_id))
"${MAIN_BIN}" -i steady_ns.yml -f \
    "main/mode=sample_generation:main/use_rom=false:sample_generation/file_path/prefix=snapshots/batch-${batch_id}/steady_ns"
