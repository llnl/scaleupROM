# Steady Navier--Stokes paper workflow

This reproduces the paper's representative 16x16, Re=25 tensorial-ROM case. It uses five component types, 2,000 random 2x2 training arrays (five jobs of 400), and 40 velocity, pressure, and supremizer modes per component. Meshes and HDF5 configurations are generated at run time and are not stored in the repository.

Install and build on Dane, then run from the configured build example:

```bash
./script/install_dane.bash
cd build/examples/steady_ns
./submit_workflow.sh
```

The submission chain runs mesh/configuration generation, sampling, `train_rom` POD and supremizer training, tensorial operator construction, a 16x16 FOM solve, and the dependent ROM comparison. If `train_rom` writes the POD bases but reaches its wall-time limit during supremizer enrichment, resume that step with `sbatch jobs/03_auxiliary.sbatch`; this fallback job is not part of the normal chain. Jobs use the paper/legacy Dane settings (`pbatch`, account `sosu`; one 2 GB single-thread task for online runs; 5 nodes and 36 MPI ranks for POD training). Override Slurm options when needed, for example `./submit_workflow.sh --account=my_account`.

Outputs are under `snapshots/`, `basis/`, `rom/`, and `results/`; logs are under `logs/`. After completion, `python3 report_results.py` prints FOM/ROM times and velocity/pressure relative errors. Each `jobs/*.sbatch` file may also be submitted separately.
