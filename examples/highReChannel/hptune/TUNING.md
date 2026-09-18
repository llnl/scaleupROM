# Tuning specification: corner2d high-Re channel

## Goal

Reach a converged corner2d channel-flow run at `nu <= 2.0e-3` that integrates at least 21.0
time units, about two flow-through times. All three criteria must hold in one trial.

`simulation_time` is that trial's own elapsed time,
`(number_of_timesteps - restart_timestep) * timestep_size` — not the total accumulated across
the restart chain.

## Notes

- **You will not get there from the start.** At the target `nu` the initial flow is under-resolved on
  this mesh and diverges. Start with a large `nu` and few timesteps. When a trial succeeds,
  gradually decrease `nu` and increase `number_of_timesteps`, restarting from a checkpoint a
  previous trial wrote.
- `number_of_timesteps` is the final timestep index, not a count, so a trial only integrates
  when `restart_timestep < number_of_timesteps`.
- `restart_timestep` must name a checkpoint that exists. They are written at step 0 and every
  `save_solution/restart_interval` (500) steps, which is why both integer parameters are on a
  step-500 grid. `restart_timestep: 0` starts from rest.
- `timestep_size` should satisfy CFL < 1.0, reported in the job log, but the actual stable
  timestep size may have a stricter condition than that.
- Wall time of 1h or less → `pdebug`; anything longer → `pbatch`.
- Cost scales with `(number_of_timesteps - restart_timestep)`. Time a short run and extrapolate
  before bidding for a long one.

## Reading list

`decide_next_trial` gets the job log as well: the CFL trend across a run is the cheapest read
on whether a transient is washing off after a change in `nu`.

The source tree is listed so a session can check what a config key does before proposing it —
`UnsteadyNSSolver::SetupInitialCondition` for the restart semantics and `SanityCheck` for what
counts as a crash.

## Specification

```yaml
goal:
  primary: nu
  criteria:
    # converged is written as H5T_NATIVE_HBOOL, which reads back as uint8 1/0 rather than a
    # bool, so this compares numerically instead of with {is: true}.
    converged:       {min: 1}
    nu:              {max: 2.0e-3}
    simulation_time: {min: 21.0}

budget:
  max_trials: 100
  max_nodes: 1
  # The job is submitted for run + decision_reserve (30m), so 23h is the largest bid that
  # still fits a 24-hour partition limit.
  max_walltime: 23h

parameters:
  solver.restart_timestep:
    int_range: [0, 99999999]
    step: 500
  time-integration.number_of_timesteps:
    int_range: [0, 99999999]
    step: 500
  time-integration.timestep_size:
    range: [1.0e-4, 3.0e-3]
  single_run.channel_flow.nu:
    range: [1.0e-3, 2.0e+2]
    scale: log

  resources.nodes:
    fixed: 1
  resources.walltime:
    choice: [30m, 1h, 2h, 4h, 8h, 16h, 23h]
  resources.queue:
    choice: [pdebug, pbatch]

allow_list:
  # PLACEHOLDER: replace /ABS/PATH/TO/scaleupROM with the real absolute path.
  interpret_run:
    - trials/{trial}/job.*.out
    - /ABS/PATH/TO/scaleupROM
  decide_next_trial:
    - trials/{trial}/job.*.out
    - /ABS/PATH/TO/scaleupROM
```
