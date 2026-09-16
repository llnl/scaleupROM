#!/bin/bash
# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT
#
# Drives the long high-Re channel runs (100,000 timesteps each).
#
# Note: assumes this is run from the build/examples/highReChannel/ directory,
# after setup_highrechannel.sh has produced the .msh.mfem meshes.
#
#   ./run_production.sh Re500               # fresh run, 1 rank
#   ./run_production.sh Re1e5 -n 8          # 8 MPI ranks
#   ./run_production.sh Re500 --resume      # continue from the newest restart file
#   ./run_production.sh Re1e5 --dt 2.0e-3   # try the nominal CFL limit

set -euo pipefail

usage() {
    cat <<'EOF'
usage: run_production.sh <Re500|Re1e5> [options]

  -n, --np N       number of MPI ranks (default: 1)
      --resume     restart from the newest restart/bump2d_<case>_*.h5
      --dt X       override time-integration/timestep_size
      --steps N    override time-integration/number_of_timesteps
      --dry-run    print the command instead of running it
EOF
}

CASE=""
NP=1
RESUME=0
DT=""
STEPS=""
DRY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        Re500|Re1e5)  CASE="$1"; shift ;;
        -n|--np)      NP="$2"; shift 2 ;;
        --resume)     RESUME=1; shift ;;
        --dt)         DT="$2"; shift 2 ;;
        --steps)      STEPS="$2"; shift 2 ;;
        --dry-run)    DRY=1; shift ;;
        -h|--help)    usage; exit 0 ;;
        *)            echo "unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z ${CASE} ]]; then
    echo "error: no case given" >&2
    usage >&2
    exit 1
fi

MAIN=../../bin/main
CONFIG="bump2d.${CASE}.yml"

# --- preflight -------------------------------------------------------------

[[ -x ${MAIN} ]] || { echo "error: ${MAIN} not found; build scaleupROM first."; exit 1; }
[[ -f ${CONFIG} ]] || { echo "error: ${CONFIG} not found."; exit 1; }

MESH=$(sed -n 's/^[[:space:]]*filename:[[:space:]]*//p' "${CONFIG}" | head -1)
if [[ ! -f ${MESH} ]]; then
    echo "error: ${CONFIG} references ${MESH}, which does not exist."
    echo "Run ./setup_highrechannel.sh first."
    exit 1
fi

mkdir -p paraview restart logs

# --- forced-input overrides ------------------------------------------------
#
# bin/main takes -f as colon-separated key=value pairs, and InputParser requires
# exactly one '=' per pair (src/input_parser.cpp). So no value here may contain
# ':' or '=' -- keep that in mind before adding a path-valued option.

FORCED=()

if [[ ${RESUME} -eq 1 ]]; then
    # Restart files are written as "%s/%s_%08d.h5", zero-padded, so the last
    # entry in lexical order is the newest step.
    LATEST=$(ls restart/bump2d_"${CASE}"_*.h5 2>/dev/null | sort | tail -1 || true)
    [[ -n ${LATEST} ]] || { echo "error: no restart file for ${CASE} in restart/."; exit 1; }
    echo "resuming from ${LATEST}"
    FORCED+=("solver/use_restart=true" "solver/restart_file=${LATEST}")
fi

[[ -n ${DT} ]] && FORCED+=("time-integration/timestep_size=${DT}")
[[ -n ${STEPS} ]] && FORCED+=("time-integration/number_of_timesteps=${STEPS}")

CMD=(mpirun -np "${NP}" "${MAIN}" -i "${CONFIG}")
if [[ ${#FORCED[@]} -gt 0 ]]; then
    CMD+=(-f "$(IFS=:; echo "${FORCED[*]}")")
fi

# --- launch ----------------------------------------------------------------

if [[ ${DRY} -eq 1 ]]; then
    printf '%q ' "${CMD[@]}"; echo
    exit 0
fi

LOG="logs/${CASE}.$(date +%Y%m%d-%H%M%S).log"
echo "running: ${CMD[*]}"
echo "logging to ${LOG}"
echo

SECONDS=0
set +e
"${CMD[@]}" 2>&1 | tee "${LOG}"
STATUS=${PIPESTATUS[0]}
set -e

printf '\n%s finished with status %d after %02d:%02d:%02d\n' \
    "${CASE}" "${STATUS}" $((SECONDS / 3600)) $((SECONDS % 3600 / 60)) $((SECONDS % 60)) \
    | tee -a "${LOG}"

exit "${STATUS}"
