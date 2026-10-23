#!/usr/bin/env bash
# ROM scaling benchmark. Run from an allocation in this directory.
set -euo pipefail

cd "$(dirname "$0")"

exe=${EXE:-../../bin/main}
deck=${DECK:-benchmark.yml}
repeats=${REPEATS:-3}
ranks=${RANKS:-"1 2 4 8"}
results=${RESULTS:-benchmarks/results.csv}
launcher=${LAUNCHER:-srun}
basis_size=${BASIS_SIZE:-10}

usage() {
  echo "Usage:"
  echo "  $0 generate SIZE...                 # SIZE is NxM or a block count"
  echo "  $0 throughput SIZE...               # one rank, varying block count"
  echo "  $0 strong SIZE [RANKS...]           # fixed problem, varying ranks"
  echo "  $0 weak BLOCKS_PER_RANK [RANKS...]  # fixed blocks per rank"
  echo "Examples:"
  echo "  $0 generate 32x32 64x64 128x128"
  echo "  $0 throughput 32x32 64x64 128x128"
  echo "  $0 strong 128x128 1 2 4 8 16"
  echo "  $0 weak 1024 1 2 4 8 16"
}

[[ $# -gt 0 ]] || { usage; exit 2; }
mode=$1
shift

if [[ $mode == generate ]]; then
  [[ $# -gt 0 ]] || { usage; exit 2; }
  exec python3 generate_benchmark.py "$@"
fi

[[ -x $exe ]] || { echo "Executable not found: $exe" >&2; exit 1; }
operator=poisson_comp
if [[ ! -f poisson_comp.h5 && -f configs/samples/poisson_comp.h5 ]]; then
  operator=configs/samples/poisson_comp
fi
[[ -f ${operator}.h5 ]] || { echo "ROM operator poisson_comp.h5 not found" >&2; exit 1; }

mkdir -p "$(dirname "$results")" benchmarks/logs
if [[ ! -f $results ]]; then
  echo "mode,nx,ny,blocks,rom_dofs,ranks,repeat,assemble_s,solve_s,kernel_s,wall_s,kernel_blocks_per_s,wall_blocks_per_s,kernel_rom_dofs_per_s,wall_rom_dofs_per_s" > "$results"
fi

run_case() {
  local kind=$1 nx=$2 ny=$3 np=$4 rep=$5
  local blocks=$((nx * ny)) rom_dofs=$((nx * ny * basis_size))
  local config="configs/benchmarks/box-channel.${nx}x${ny}.h5"
  local log="benchmarks/logs/${kind}.${nx}x${ny}.n${np}.r${rep}.log"
  [[ -f $config ]] || python3 generate_structured_benchmark.py "${nx}x${ny}"

  local overrides="mesh/component-wise/global_config=${config}"
  overrides+=":visualization/enabled=false"
  overrides+=":model_reduction/compare_solution/enabled=false"
  overrides+=":model_reduction/save_linear_system/enabled=false"
  overrides+=":model_reduction/save_operator/prefix=${operator}"

  echo "[$kind] ${nx}x${ny} (${blocks} blocks), ${np} rank(s), repeat ${rep}/${repeats}"
  local start_ns end_ns wall
  start_ns=$(date +%s%N)
  OMP_NUM_THREADS=1 "$launcher" -n "$np" --cpu-bind=cores \
    "$exe" -i "$deck" -f "$overrides" > "$log" 2>&1
  end_ns=$(date +%s%N)
  wall=$(awk -v a="$start_ns" -v b="$end_ns" 'BEGIN {printf "%.9g", (b-a)/1e9}')

  # With multiple ranks, use the slowest reported phase.
  local assemble solve total rate wall_rate dof_rate wall_dof_rate
  assemble=$(awk '/ROM-assemble time:/ {if ($3+0>m) m=$3} END {print m+0}' "$log")
  solve=$(awk '/ROM-solve time:/ {if ($3+0>m) m=$3} END {print m+0}' "$log")
  total=$(awk -v a="$assemble" -v s="$solve" 'BEGIN {printf "%.9g", a+s}')
  rate=$(awk -v n="$blocks" -v t="$total" 'BEGIN {printf "%.9g", n/t}')
  wall_rate=$(awk -v n="$blocks" -v t="$wall" 'BEGIN {printf "%.9g", n/t}')
  dof_rate=$(awk -v n="$rom_dofs" -v t="$total" 'BEGIN {printf "%.9g", n/t}')
  wall_dof_rate=$(awk -v n="$rom_dofs" -v t="$wall" 'BEGIN {printf "%.9g", n/t}')
  echo "$kind,$nx,$ny,$blocks,$rom_dofs,$np,$rep,$assemble,$solve,$total,$wall,$rate,$wall_rate,$dof_rate,$wall_dof_rate" | tee -a "$results"
}

split_size() {
  local spec=$1 blocks nx ny
  if [[ $spec == *x* ]]; then
    nx=${spec%x*}; ny=${spec#*x}
  else
    blocks=$spec; nx=$(awk -v n="$blocks" 'BEGIN {printf "%d", sqrt(n)}')
    while (( blocks % nx != 0 )); do ((nx--)); done
    ny=$((blocks / nx))
  fi
  echo "$nx $ny"
}

case $mode in
  throughput)
    [[ $# -gt 0 ]] || { usage; exit 2; }
    for size in "$@"; do read -r nx ny <<< "$(split_size "$size")"; for ((r=1; r<=repeats; r++)); do run_case throughput "$nx" "$ny" 1 "$r"; done; done
    ;;
  strong)
    [[ $# -gt 0 ]] || { usage; exit 2; }
    read -r nx ny <<< "$(split_size "$1")"; shift
    test_ranks=${*:-$ranks}
    for np in $test_ranks; do for ((r=1; r<=repeats; r++)); do run_case strong "$nx" "$ny" "$np" "$r"; done; done
    ;;
  weak)
    [[ $# -gt 0 ]] || { usage; exit 2; }
    per_rank=$1; shift
    test_ranks=${*:-$ranks}
    for np in $test_ranks; do read -r nx ny <<< "$(split_size "$((per_rank * np))")"; for ((r=1; r<=repeats; r++)); do run_case weak "$nx" "$ny" "$np" "$r"; done; done
    ;;
  *) usage; exit 2 ;;
esac

echo "Results appended to $results"
