# High-Reynolds channel flow over a bump

A monolithic full-order `UnsteadyNSSolver` reproduction of the OpenLB `bump2d`
applications in [`openlb_examples/`](../../openlb_examples), at the two Reynolds
numbers those runs use: **Re = 500** (`bump2d_laminar.cpp`) and **Re = 1e5**
(`bump2d_turbulent.cpp`).

OpenLB is a lattice Boltzmann code, so nothing about its discretization carries
over directly. Most of this file is the reasoning behind *how* the two setups were
made comparable — that is the part worth checking before trusting a comparison.

```
   y=1.5 +-------------------------------------------------+   1  inflow      (x = 0)
         |                                                 |   2  no-slip wall
   y=.75 +-----------+###########+-------------------------+   3  outflow     (x = 25)
         |           |###########|                         |
   y=0   +-----------+###########+-------------------------+
        x=0         x=5        x=10                      x=25
```

Fluid domain `[0, 25] x [0, 1.5]` minus the step `[5, 10] x [0, 0.75]`, identical
to the OpenLB geometry. Inflow is a Poiseuille profile of peak 1.5 (mean 1.0 =
OpenLB's `charPhysVelocity`); the outflow is traction-free.

## Quick start

From `build/examples/highReChannel/`:

```bash
./setup_highrechannel.sh          # generate meshes, convert to MFEM, make output dirs
../../bin/main -i bump2d.smoke.yml   # 200-step smoke test, under a minute
./run_production.sh Re500 -n 8    # the real runs
./run_production.sh Re1e5 -n 8
```

`setup_highrechannel.sh` is a prerequisite, not an optional regeneration step: the
meshes are not committed (`*.msh` and `*.msh.mfem` are in the repo's
`.gitignore`), so the configs have nothing to load until it has run once. It needs
the gmsh **Python module** (`pip3 install gmsh`); note that the `gmsh` apt package
installs the CLI only, which is enough for CMake's `find_program(GMSH gmsh)` — and
hence for `utils/gmsh2mfem` to be built — but not for `generate_mesh.py`.

## The mapping

| | OpenLB laminar | OpenLB turbulent | `bump2d.Re500` | `bump2d.Re1e5` |
|---|---|---|---|---|
| domain | `[0,25]x[0,1.5]` less step | identical | identical | identical |
| viscosity `nu` | 2e-3 | 1e-5 | 2e-3 | 1e-5 |
| Reynolds number | 500 | 1e5 | 500 | 1e5 |
| resolution | `dx = 1/60` | `dx = 1/30` | `h = 1/32`, `h_eff = 1/64` | `h = 1/16`, `h_eff = 1/32` |
| cells / elements | 121,500 | 30,375 | 34,560 quads | 8,640 quads |
| unknowns | 9 `f_i` per cell | 9 `f_i` per cell | ~760k DOF | ~190k DOF |
| timestep | 3.33e-4 | 6.67e-4 | 1.0e-3 | 1.0e-3 |
| advective CFL | 0.02 | 0.02 | 0.19 | 0.10 |
| simulated time | 100 | 100 | 100 | 100 |
| turbulence model | none | Smagorinsky, Cs=0.15 | none | none (implicit, via LF flux) |
| startup | rest + 20s ramp | rest + 20s ramp | consistent Poiseuille IC | consistent Poiseuille IC |

### Resolution: match node spacing, not element size

D2Q9 carries three macroscopic unknowns per cell on a lattice of spacing `dx`. The
invariant worth matching is the spacing between *velocity nodes*, which is also
what the solver's own CFL diagnostic uses — `UnsteadyNSSolver::ComputeCFL` divides
the element size by the element order. So:

```
h_eff = h / p_u,    p_u = uorder = discretization/order + 1
```

With `discretization/order: 1` the velocity space is Q2, so `p_u = 2` and
`h = 2 dx`. That gives `h = 1/32` to match `dx = 1/60`, and `h = 1/16` to match
`dx = 1/30`.

The comparison is not being quietly tilted: DG Q2/Q1 stores `2*9 + 4 = 22` DOF per
element of area `(2 dx)^2`, i.e. **5.5 DOF per `dx^2`**, against D2Q9's **9
distributions per `dx^2`**.

### Why the higher-Re case gets the *coarser* mesh

This looks backwards and is deliberate — it mirrors OpenLB, where the two source
files have byte-identical geometry and differ only in `N` (60 vs 30). OpenLB
coarsens the harder problem for three reasons:

1. **Cost compounds.** `dt = 1/(M*N)`, so `N` sets both the grid and the step
   count. Going 60 -> 30 is 4x fewer cells *and* 2x fewer steps: 8x cheaper for the
   same 100 units of simulated time.
2. **Refining buys no physics.** At Re = 1e5 the first cell centre sits at
   `y+ ~ 150` at N=30 and `y+ ~ 75` at N=60. Wall-resolved LES wants `y+ <~ 1`, so
   both grids are equally on the wrong side of that threshold — 2x refinement
   changes cost, not regime.
3. **Coarser is more stable, with an SGS model.** Lattice viscosity is
   `nu_LB = nu * N/M` = 6e-6 at N=30, giving a relaxation time `tau = 0.500018`:
   molecular viscosity contributes nothing, and the Smagorinsky term is doing all
   the work. Since `nu_t,LB = Cs^2 |S_LB|` scales as `1/(M*N)`, the coarser lattice
   carries *more* subgrid dissipation and sits further from the unstable
   `tau = 1/2` limit. Refining the turbulent run would be more expensive **and**
   less stable.

The direct consequence: `bump2d.Re500` is the memory-heavy run here (~760k DOF),
not `bump2d.Re1e5`.

### Timestep: set by explicit convection

`UnsteadyNSSolver::Step` is semi-implicit. The BDF1 mass term, the viscous
operator and the pressure/divergence coupling are all implicit, inside the
MUMPS-factorized saddle-point system. Convection is explicit (AB1 / forward
Euler). So:

- there is **no viscous timestep limit** — it would have been `h_eff^2/(2 nu)`,
  not binding at either Reynolds number anyway;
- there is **no Mach limit**, which is the whole point of not being an LBM.
  OpenLB's `dt` is pinned by `u_LB = charU * dt/dx = 1/M = 0.02`, independent of
  `nu` — it is a stability constraint, not an accuracy choice, and there is
  nothing in it to "match";
- what remains is the advective CFL on the explicit term. Forward Euler with a
  central operator is unconditionally unstable; this is stable only because
  `DGLaxFriedrichsFluxIntegrator` supplies `|u.n|`-scaled upwind dissipation.
  Budget `CFL ~ 0.2`.

The channel narrows 1.5 -> 0.75 over the step, so by continuity the mean velocity
doubles and `u_max ~ 3.0` in the gap. Then
`dt = CFL * h_eff / u_max`:

- **Re = 500**: `0.2 * (1/64) / 3 = 1.04e-3` -> `dt = 1.0e-3`, CFL ~ 0.19.
- **Re = 1e5**: the same formula gives 2.1e-3 on the coarser mesh, but this config
  uses `dt = 1.0e-3` (CFL ~ 0.10), halved for margin against steeper near-wall
  gradients and chaotic excursions above `u_max = 3`.

Both configs therefore run at `dt = 1.0e-3` — for different reasons, one
CFL-limited and one margin-limited, not by copy-paste. If the Re=1e5 case proves
stable, raising it to 2.0e-3 halves its wall time and is the first lever to pull.

For scale: our advective CFL is ~0.19 against OpenLB's 0.02, so at matched spatial
resolution the timestep is only about 3x larger. Landing in the same ballpark is
the sanity check that the mapping is right — it is not a speedup claim.

### What is *not* matched

**No turbulence model at Re = 1e5.** scaleupROM has no SGS model, so the Re=1e5
run is an under-resolved calculation stabilized by numerical dissipation: the
interior Lax–Friedrichs flux plays the role OpenLB gives to Smagorinsky. That is
an honest counterpart — the OpenLB run's molecular viscosity is negligible too
(`tau = 0.500018`) — but it is an implicit-LES comparison, so expect qualitative
agreement in separation and shedding behaviour, not pointwise agreement. In 2D
there is no vortex stretching and no real turbulent cascade on either side.

This is why both configs set `full-discrete-galerkin: true`: the interior
Lax–Friedrichs integrator is only added when `full_dg` is on
(`SteadyNSSolver::BuildDomainOperators`). Without it the interior convection
operator is unstabilized central, which is fatal at a cell Reynolds number of
~3e3. At Re=500 the cell Reynolds number is only ~8, so that config may safely
fall back to continuous Taylor–Hood (`full-discrete-galerkin: false`, ~310k DOF)
if the MUMPS factorization does not fit in memory. **Do not do this at Re = 1e5** —
it removes the only stabilization in the scheme.

**No inflow ramp.** OpenLB starts from rest and ramps the inlet over the first 20
units via `PolynomialStartScale`. The `backward_facing_step` problem's initial
condition is instead the Poiseuille profile itself, which is divergence-free and
*consistent with the inflow BC at t = 0*, so no ramp is needed. The only initial
inconsistency is no-slip on the step faces, which relaxes within a few steps as a
classic impulsive start. Comparable windows are therefore OpenLB `t in [20, 100]`
against ours `t in [0, 80]`.

## Configuration notes

The problem is [`backward_facing_step`](../../src/parameterized_problem.cpp), used
unmodified. Its three boundary attributes map onto the mesh's physical curve tags:
1 = Dirichlet inflow, 2 = zero (no-slip wall), 3 = Neumann outflow — the same
convention as [`utils/gmsh/bfs.geo`](../../utils/gmsh/bfs.geo).

- `u0` is the Poiseuille **peak**, not the mean. `u0: 1.5` gives mean 1.0.
- `y0`/`y1` bound the inflow profile; here the inlet spans the full channel, so
  `0.0` and `1.5`.
- The viscosity that takes effect is `single_run/backward_facing_step/nu`, which
  overwrites `stokes/nu` during `SetParameterizedProblem`. Both are set to the
  same value in these configs to avoid confusion.
- `amp0`/`ky0`/`freq0`/`t_offset0` add a time-periodic transverse perturbation at
  the inlet. Off by default; try `amp0: 0.01, ky0: 1.0, freq0: 0.5` to trip the
  flow if it stays stubbornly steady.
- The Neumann outflow means `pres_dbc` is true, so the pressure-constant removal
  and complementary-flux paths stay off — matching OpenLB's pressure outflow BC.

## Meshes

`generate_mesh.py` builds the domain from five rectangular blocks that share whole
curves, so transfinite meshing is conforming by construction and **any** element
size works — there is no divisibility constraint:

```bash
python3 generate_mesh.py --h 0.0625 --check    # 8,640 quads
python3 generate_mesh.py --h 0.03125 --check   # 34,560 quads
```

`--check` re-reads the written `.msh` and asserts the invariants that matter
downstream: pure quads (mixed element types trip a `SteadyNSSolver` assertion),
conformity, no orphan nodes, exactly one surface tag (so the mesh loads as a
single subdomain), exactly three boundary tags on the correct geometry, and that
the step is genuinely a hole.

The script also prints the realized `dx`, `h_eff`, and the suggested
`timestep_size` for CFL = 0.2, so the numbers in this README can be re-derived for
any resolution.

## If the smoke test disagrees

The load-bearing assumption in every `timestep_size` here is `u_max ~ 3.0`. The
smoke test reports CFL every 20 steps; it should settle near 0.19. If it does not,
rescale the production `timestep_size` values by the observed ratio before
launching the long runs.

If the Re=1e5 run goes unstable, escalate in this order: halve `--dt`; refine to
`--h 0.03125`; add an inlet trip via `amp0`/`ky0`/`freq0`.
