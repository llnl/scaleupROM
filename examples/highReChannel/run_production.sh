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
#   ./run_production.sh bump2d.Re500                      # fresh run, 1 rank
#   ./run_production.sh corner2d.Re1e5 -n 8               # 8 MPI ranks
#   ./run_production.sh bump2d.Re500 --resume             # continue from newest restart
#   ./run_production.sh corner2d.Re500 --problem channel_flow   # start from rest
#   ./run_production.sh bump2d.Re1e5 --dt 2.0e-3          # try the nominal CFL limit

set -euo pipefail

usage() {
    cat <<'EOF'
usage: run_production.sh <case> [options]

  <case> is the basename of a config in this directory, e.g.
      bump2d.Re500   bump2d.Re1e5   bump2d.smoke
      corner2d.Re500 corner2d.Re1e5 corner2d.smoke

  -n, --np N          number of MPI ranks (default: 1)
      --resume        restart from the newest checkpoint for this case
      --problem NAME  backward_facing_step | channel_flow
                      (also swaps the mesh to the matching tag variant)
      --dt X          override time-integration/timestep_size
      --steps N       override time-integration/number_of_timesteps
      --dry-run       print the command instead of running it
EOF
}

CASE=""
NP=1
RESUME=0
PROBLEM=""
DT=""
STEPS=""
DRY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--np)      NP="$2"; shift 2 ;;
        --resume)     RESUME=1; shift ;;
        --problem)    PROBLEM="$2"; shift 2 ;;
        --dt)         DT="$2"; shift 2 ;;
        --steps)      STEPS="$2"; shift 2 ;;
        --dry-run)    DRY=1; shift ;;
        -h|--help)    usage; exit 0 ;;
        -*)           echo "unknown option: $1" >&2; usage >&2; exit 1 ;;
        *)            CASE="$1"; shift ;;
    esac
done

if [[ -z ${CASE} ]]; then
    echo "error: no case given" >&2
    usage >&2
    exit 1
fi

MAIN=../../bin/main
CONFIG="${CASE}.yml"

# --- config reader ---------------------------------------------------------
#
# PyYAML is guaranteed in the container (docker/Dockerfile pip-installs it), and
# is a lot less fragile than sed against nested keys that repeat, such as the
# two file_path/prefix entries.

yget() {
    python3 -c "
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
for k in sys.argv[2].split('/'):
    if not isinstance(d, dict) or k not in d:
        sys.exit(0)
    d = d[k]
print(d)
" "${CONFIG}" "$1"
}

# --- preflight -------------------------------------------------------------

[[ -x ${MAIN} ]] || { echo "error: ${MAIN} not found; build scaleupROM first."; exit 1; }
python3 -c "import yaml" 2>/dev/null || {
    echo "error: python3 with PyYAML is required to read the config."
    echo "Install it with 'pip3 install PyYAML'."
    exit 1
}
if [[ ! -f ${CONFIG} ]]; then
    echo "error: ${CONFIG} not found. Available cases:"
    for f in *.yml; do echo "    ${f%.yml}"; done
    exit 1
fi

MESH=$(yget mesh/filename)

# --- forced-input overrides ------------------------------------------------
#
# bin/main takes -f as colon-separated key=value pairs, and InputParser requires
# exactly one '=' per pair (src/input_parser.cpp). So no value here may contain
# ':' or '=' -- keep that in mind before adding a path-valued option.

FORCED=()

if [[ -n ${PROBLEM} ]]; then
    # backward_facing_step wants tags 1=in/2=wall/3=out; channel_flow wants
    # 4=in/2=out/1=wall. Same geometry, different physical group numbers, so the
    # problem and the mesh have to move together.
    case "${PROBLEM}" in
        backward_facing_step) want=bfs ;;
        channel_flow)         want=chan ;;
        *) echo "error: --problem must be backward_facing_step or channel_flow" >&2; exit 1 ;;
    esac
    MESH="${MESH%.bfs.msh.mfem}"; MESH="${MESH%.chan.msh.mfem}"
    MESH="${MESH}.${want}.msh.mfem"
    FORCED+=("parameterized_problem/name=${PROBLEM}" "mesh/filename=${MESH}")
    echo "problem: ${PROBLEM} (mesh ${MESH})"
fi

if [[ ! -f ${MESH} ]]; then
    echo "error: ${CONFIG} needs ${MESH}, which does not exist."
    echo "Run ./setup_highrechannel.sh first."
    exit 1
fi

mkdir -p paraview restart logs

if [[ ${RESUME} -eq 1 ]]; then
    SOLDIR=$(yget save_solution/file_path/directory)
    SOLPRE=$(yget save_solution/file_path/prefix)
    if [[ -z ${SOLDIR} || -z ${SOLPRE} ]]; then
        echo "error: ${CONFIG} does not enable save_solution, so there is nothing to resume from."
        exit 1
    fi
    # Restart files are written as "%s/%s_%08d.h5", zero-padded, so the last
    # entry in lexical order is the newest step.
    LATEST=$(ls "${SOLDIR}/${SOLPRE}"_*.h5 2>/dev/null | sort | tail -1 || true)
    [[ -n ${LATEST} ]] || { echo "error: no restart file matching ${SOLDIR}/${SOLPRE}_*.h5"; exit 1; }
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
