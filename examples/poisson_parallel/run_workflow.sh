#!/usr/bin/env bash
# Run the component Poisson ROM workflow from training through scaling plots.
# This script is copied to build/examples/poisson_parallel by CMake and should
# be run there, preferably inside an allocation for the benchmark stages.
set -euo pipefail

cd "$(dirname "$0")"

exe=${EXE:-../../bin/main}
samples=${SAMPLES:-50}
launcher=${LAUNCHER:-srun}
repeats=${REPEATS:-3}
ranks=${RANKS:-"1 2 4 8 16"}
results=${RESULTS:-benchmarks/results.csv}

usage() {
  cat <<'EOF'
Usage: ./run_workflow.sh [prepare|train|benchmark|plots|all]

  prepare    Generate sample configurations and create output directories.
  train      Prepare, generate FOM samples, train bases, and build the ROM.
  benchmark  Run the published throughput, strong-, and weak-scaling sweeps.
  plots      Regenerate plots from the benchmark CSV.
  all        Run every stage (default).

Environment overrides:
  EXE, SAMPLES, LAUNCHER, REPEATS, RANKS, RESULTS, MPLCONFIGDIR
EOF
}

prepare() {
  mkdir -p snapshots basis rom-system paraview benchmarks/logs benchmarks/plots
  if [[ ! -f configs/box-channel.comp.h5 || \
        ! -f configs/samples/box-channel.config-00624.h5 ]]; then
    ./setup_poisson.sh
  fi
}

train() {
  prepare
  [[ -x $exe ]] || { echo "Executable not found: $exe" >&2; exit 1; }

  echo "Generating $samples FOM training samples"
  "$exe" -i poisson.sampling.yml \
    -f "sample_generation/random_sample_generator/number_of_samples=${samples}"

  echo "Training component reduced bases"
  "$exe" -i poisson.sampling.yml \
    -f "main/mode=train_rom:sample_generation/random_sample_generator/number_of_samples=${samples}"

  echo "Building the component ROM operator"
  "$exe" -i poisson.sampling.yml -f "main/mode=build_rom:main/use_rom=true"
}

benchmark() {
  [[ -f poisson_comp.h5 || -f configs/samples/poisson_comp.h5 ]] || {
    echo "ROM operator not found; run '$0 train' first" >&2
    exit 1
  }

  EXE=$exe LAUNCHER=$launcher REPEATS=$repeats RESULTS=$results \
    ./benchmark.sh throughput 16x16 32x32 64x64 80x80 92x92 128x128
  EXE=$exe LAUNCHER=$launcher REPEATS=$repeats RESULTS=$results \
    ./benchmark.sh strong 128x128 $ranks
  EXE=$exe LAUNCHER=$launcher REPEATS=$repeats RESULTS=$results \
    ./benchmark.sh weak 80 $ranks
}

plots() {
  MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/matplotlib-cache} \
    ./plot_benchmarks.py --input "$results"
}

stage=${1:-all}
case $stage in
  prepare) prepare ;;
  train) train ;;
  benchmark) benchmark ;;
  plots) plots ;;
  all) train; benchmark; plots ;;
  -h|--help|help) usage ;;
  *) usage >&2; exit 2 ;;
esac
