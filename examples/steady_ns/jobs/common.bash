set -euo pipefail

EXAMPLE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "${EXAMPLE_DIR}"

if [[ -n ${SCALEUPROM_MAIN:-} ]]; then
    MAIN_BIN=${SCALEUPROM_MAIN}
elif [[ -x ../../bin/main ]]; then
    MAIN_BIN=../../bin/main
elif [[ -x ../../build/bin/main ]]; then
    MAIN_BIN=../../build/bin/main
else
    echo "Cannot find the scaleupROM main executable; set SCALEUPROM_MAIN." >&2
    exit 1
fi

export OMP_NUM_THREADS=1
mkdir -p logs results rom basis snapshots

