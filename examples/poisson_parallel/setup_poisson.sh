#!/bin/bash
# Generates input meshes for the Poisson example
# Note: assumes this is run from the build/examples/poisson_parallel directory
comp_script="../../../utils/python/box_comp_config.py"
if [[ ! -f ${comp_script} ]]; then
    echo "Could not find box_comp_config.py script"
    exit 1
fi

# Generate box-channel.comp.h5
python3 ${comp_script}

# Generate all sample config meshes
python3 generate_configs.py

mkdir -p configs/
mv box-channel.comp.h5 configs/
mv test.box-channel.8x8.h5 configs/

mkdir -p configs/samples/
mv *.h5 configs/samples/
