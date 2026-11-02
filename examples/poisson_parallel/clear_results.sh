#!/usr/bin/env bash
set -euo pipefail

name=${1:-old}
archive="benchmarks/plots/$name"

if [[ -e $archive ]]; then
  echo "Archive already exists: $archive" >&2
  echo "Choose another name: ./clear_results.sh NAME" >&2
  exit 1
fi

mkdir -p "$archive" benchmarks/logs
shopt -s nullglob
plots=(benchmarks/plots/*.png)
if ((${#plots[@]})); then
  mv "${plots[@]}" "$archive/"
fi

rm -f benchmarks/results.csv
rm -f benchmarks/profile_summary.csv benchmarks/profile_exclusive_summary.csv
rm -f benchmarks/logs/*.log

echo "Results cleared; old plots saved in $archive"
