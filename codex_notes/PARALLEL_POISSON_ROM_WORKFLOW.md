# Parallel Poisson ROM Reproduction Notes

Run these commands from:

```text
build/examples/poisson_parallel
```

Commands are shown on one line. Serial steps invoke the executable directly;
parallel steps will use `srun` rather than `mpirun`.

The complete workflow is automated by:

```bash
./run_workflow.sh all
```

Run it inside an allocation. Individual stages are available as `prepare`,
`train`, `benchmark`, and `plots`; run `./run_workflow.sh --help` for environment
overrides. The commands below document what the training stage does.

## 1. Create output directories

```bash
mkdir -p snapshots basis rom-system paraview/8x8
```

Creates the directories referenced by the example inputs. The supplied setup
script does not create them. Snapshot generation uses `snapshots`; later stages
use `basis`, `rom-system`, and `paraview/8x8`.

## 2. Generate 50 FOM training samples

```bash
../../bin/main -i poisson.sampling.yml -f "sample_generation/random_sample_generator/number_of_samples=50"
```

Runs 50 serial Poisson FOM samples using randomly selected 2-by-2 component
layouts and forcing/boundary parameters. Component solution snapshots are grouped
by component type and written to `snapshots/`. This generates training data; it
does not train or run a ROM.

## 3. Train component reduced bases

```bash
../../bin/main -i poisson.sampling.yml -f "main/mode=train_rom:sample_generation/random_sample_generator/number_of_samples=50"
```

Reads the component snapshot databases, performs an SVD/POD separately for each
component type, and saves ten basis vectors per type under `basis/`. The printed
singular values measure the importance of successive solution patterns; the
energy fraction reports how much snapshot variation is captured cumulatively by
the retained patterns.


## 4. Build the component ROM operator

```bash
../../bin/main -i poisson.sampling.yml -f "main/mode=build_rom:main/use_rom=true"
```

Reads the trained component bases and writes `poisson_comp.h5`, which contains
the reduced component operators consumed by the scaling benchmarks.

## 5. Run all published scaling campaigns

```bash
./run_workflow.sh benchmark
```

This runs the single-rank throughput sweep, the fixed 128-by-128 strong-scaling
sweep, and the 80-blocks-per-rank weak-scaling sweep. Each case is repeated
three times on 1, 2, 4, 8, and 16 ranks where applicable.
