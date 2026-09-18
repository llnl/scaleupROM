#!/bin/sh
# Start the campaign with: trial_epilogue --init --study <study>

#SBATCH --job-name=hptune-corner2d
#SBATCH --account=PLACEHOLDER_BANK
# No --nodes/--ntasks/--time/--partition/--output/--export: the loop passes those.

set -u

module load mvapich2

: "${STUDY:?STUDY must be set}"
: "${VENV:=$STUDY/venv}"
if [ -d "$VENV" ]; then
    # shellcheck disable=SC1091
    . "$VENV/bin/activate"
fi

# Checked here rather than at the handoff on the last line, so a missing venv costs a queue
# slot instead of a queue slot plus the whole run.
if ! command -v trial_epilogue >/dev/null 2>&1; then
    echo "job.sbatch: hptune is not on PATH. Pass VENV=<study>/venv." >&2
    exit 127
fi

: "${TRIAL:?TRIAL must be set}"
: "${HPTUNE_RUN_TIMELIMIT:=00:30:00}"

TRIAL_DIR="$STUDY/trials/$TRIAL"
mkdir -p "$TRIAL_DIR"

# main resolves every relative path against its cwd -- mesh, restart directory, metrics file --
# so the run happens in the study directory, not the trial directory. That is deliberate: the
# restart files are shared state, and a later trial resumes from an earlier trial's checkpoint.
cd "$STUDY" || exit 1
mkdir -p restart                 # H5Fcreate does not create its own directory

# A run killed by the time limit writes no metrics file, and H5Fcreate truncates rather than
# appends, so clear any leftover before the run rather than after.
rm -f "$STUDY/metrics.h5"

STARTED=$(date +%s)

srun -N "$NNODES" -n "$NPROCS" -t "$HPTUNE_RUN_TIMELIMIT" \
    "$STUDY/main" -i "$TRIAL_DIR/config.yaml"
EXIT_STATUS=$?

# main cannot be told to write the metrics file anywhere but its cwd.
if [ -f "$STUDY/metrics.h5" ]; then
    mv "$STUDY/metrics.h5" "$TRIAL_DIR/metrics.h5"
fi

exec trial_epilogue \
    --study "$STUDY" --trial "$TRIAL" \
    --exit-status "$EXIT_STATUS" --started-at "$STARTED"
