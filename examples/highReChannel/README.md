# High-Reynolds channel flows

Monolithic full-order `UnsteadyNSSolver` reproductions of two OpenLB
applications in [`openlb_examples/`](../../openlb_examples), each at the two
Reynolds numbers those runs use — **Re = 500** (`*_laminar.cpp`) and **Re = 1e5**
(`*_turbulent.cpp`):

- **`bump2d`** — a rectangular step on the floor of a long channel.
- **`corner2d`** — an elbow: a duct that expands into a square chamber, turns
  through 90 degrees, and contracts back into a duct.

OpenLB is a lattice Boltzmann code, so nothing about its discretization carries
over directly. Most of this file is the reasoning behind *how* the setups were
made comparable — that is the part worth checking before trusting a comparison.

```
bump2d:  [0,25] x [0,1.5] less the step [5,10] x [0,0.75]

   y=1.5 +-------------------------------------------------+
         |                                                 |
   y=.75 +-----------+###########+-------------------------+
         |           |###########|                         |
   y=0   +-----------+###########+-------------------------+
        x=0         x=5        x=10                      x=25

corner2d: [0,11] x [0,5] less [0,1]x[0,4] and [6,11]x[1,5]

   y=5 +---+-----------+###############+
       |   |           |###############|      inflow   x = 0
   y=4 +---+-----------+###############+      outflow  x = 11 (bump2d: 25)
       |###|           |###############|      walls    everything else
   y=1 |###+-----------+---------------+
       |###|                           |
   y=0 +---+---------------------------+
      x=0 x=1         x=6            x=11
```

Inflow is a Poiseuille profile of peak 1.5 — mean exactly 1.0, which is OpenLB's
`charPhysVelocity` — across the inlet duct; the outflow is traction-free. For
bump2d the inlet spans the full channel height `[0, 1.5]`; for corner2d it is the
`y in [4, 5]` duct.

## Quick start

From `build/examples/highReChannel/`:

```bash
./setup_highrechannel.sh                  # meshes -> MFEM, plus output dirs
../../bin/main -i bump2d.smoke.yml        # 200-step smoke test, under a minute
../../bin/main -i corner2d.smoke.yml
./run_production.sh bump2d.Re500 -n 8     # the real runs
./run_production.sh bump2d.Re1e5 -n 8
./run_production.sh corner2d.Re500 -n 8
./run_production.sh corner2d.Re1e5 -n 8
```

`setup_highrechannel.sh` is a prerequisite, not an optional regeneration step: the
meshes are not committed (`*.msh` and `*.msh.mfem` are in the repo's
`.gitignore`), so the configs have nothing to load until it has run once. It needs
the gmsh **Python module** (`pip3 install gmsh`); note that the `gmsh` apt package
installs the CLI only, which is enough for CMake's `find_program(GMSH gmsh)` — and
hence for `utils/gmsh2mfem` to be built — but not for `generate_mesh.py`.

## The mapping

Both OpenLB cases use `M = 50`, `charPhysLength = charPhysVelocity = 1`,
`maxPhysT = 100`, and `N = 60` (laminar) / `N = 30` (turbulent), so the mapping
below is the same recipe applied to two geometries.

| | OpenLB laminar | OpenLB turbulent | `*.Re500` | `*.Re1e5` |
|---|---|---|---|---|
| viscosity `nu` | 2e-3 | 1e-5 | 2e-3 | 1e-5 |
| Reynolds number | 500 | 1e5 | 500 | 1e5 |
| resolution | `dx = 1/60` | `dx = 1/30` | `h = 1/32`, `h_eff = 1/64` | `h = 1/16`, `h_eff = 1/32` |
| timestep | 3.33e-4 | 6.67e-4 | 1.0e-3 | 1.0e-3 |
| simulated time | 100 | 100 | 100 (100,000 steps) | 100 (100,000 steps) |
| turbulence model | none | Smagorinsky, Cs=0.15 | none | none (implicit, via LF flux) |
| startup | rest + 20s ramp | rest + 20s ramp | see *Switching the inflow problem* | " |

Per geometry:

| | bump2d Re500 | bump2d Re1e5 | corner2d Re500 | corner2d Re1e5 |
|---|---|---|---|---|
| LBM cells | 121,500 | 30,375 | 111,600 | 27,900 |
| element budget | 34,560 quads | 8,640 quads | 31,744 quads | 7,936 quads |
| achieved (graded) | 34,817 | 8,756 | 31,412 | 8,129 |
| DG unknowns | ~766k | ~193k | ~691k | ~179k |
| `h_min` .. `h_max` | .0164 .. .0493 | .0332 .. .0997 | .0156 .. .0469 | .0313 .. .0938 |
| `u_max` budget | 3.0 | 3.0 | 2.5 | 2.5 |
| advective CFL, uniform | 0.19 | 0.10 | 0.16 | 0.08 |
| advective CFL, graded | ~0.36 | ~0.18 | ~0.32 | ~0.16 |

`u_max` differs because bump2d's channel narrows 1.5 → 0.75 over the step, so by
continuity the mean velocity doubles and the peak reaches ~3.0 in the gap.
corner2d has no area contraction — inlet and outlet are both height 1 — so the
bulk peak stays ~1.5, and the 2.5 budget is margin for local speed-up at the two
convex corners `(1,4)` and `(6,1)`.

For reference, OpenLB's own advective CFL is `u_LB = 1/M = 0.02` in every case.

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

The meshes are **graded**, so `h` is no longer a single number: `--h` now sets the
DOF *budget* (the element count a uniform mesh of that size would have had) while
the actual spacing runs from `h_min` at walls to `3 h_min` in the core. See
*Wall-distance grading* below. Element counts are within 3% of the uniform
figures, so the DOF and cost numbers in this file still hold.

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

### Wall-distance grading

Uniform spacing is wasteful: at Re = 1e5 the boundary layer at the step is
`δ ≈ 0.134`, and uniform `h_eff = 1/32` puts only **4.3 velocity nodes across
it**, while mid-channel is over-resolved for what happens there. The meshes are
therefore graded with `h_max/h_min = 3` **at fixed element count**, so this is a
redistribution of resolution, not an increase in cost.

Size is driven by a gmsh distance-to-wall field:

```
Field[1] = Distance,  CurvesList = <wall curves only>
Field[2] = Threshold, SizeMin = h_min,  SizeMax = 3 h_min
                      DistMin = 0,      DistMax = layer
```

Because the field measures true Euclidean distance to the wall *curves*, the
step's **vertical faces** and the **convex corners** are refined exactly like the
floor and ceiling, and the inflow/outflow curves — simply absent from
`CurvesList` — attract no refinement at all. Measured on `bump2d` at `--h 0.0625`
(`h_min = 0.0332`, `h_max = 0.0997`):

| probe | distance to nearest wall | element size |
|---|---|---|
| step vertical face `(5, 0.40)` | 0.000 | 0.0327 |
| corner `(5, 0.75)` | 0.020 | 0.0387 |
| floor `(15, 0)` / ceiling `(15, 1.5)` | 0.020 | 0.0366 / 0.0375 |
| interior, mid-channel `(2.5, 0.75)` | 0.750 | **0.1115** |
| interior, core `(18, 0.375)` | 0.375 | **0.0954** |

`layer` (`DistMax`) is tied to geometry, not to `h`: half the characteristic gap
height — 0.375 for bump2d, 0.5 for corner2d — so the size reaches `h_max` exactly
at mid-gap. That makes the per-element growth rate come out at **≈ 1.21**, which
is textbook boundary-layer practice (1.05–1.2).

**Why the mesh is unstructured.** The original transfinite mesher cannot express
this. In a conforming structured mesh the node distribution along an axis is
shared by every block in the perpendicular strip, so refinement cannot be
localized to *part* of a line: the step's left face is a wall only for
`y ∈ [0, 0.75]`, but clustering toward `x = 5` would unavoidably cluster it for
`y ∈ [0.75, 1.5]` too, and likewise `y = 0.75` would be refined across the whole
channel rather than just over the step. Grading with no spurious interior
refinement would have left only `y = 0` and `y = 1.5` — a pure `f(y)` law that
ignores the vertical faces entirely.

**What grading does not fix.** At Re = 500 the mesh was already wall-resolved
(`y+₁ ≈ 0.64`), so grading there only buys a coarser core. At Re = 1e5,
`y+₁ ≈ 143` uniform and ≈ 72 graded, against the `y+₁ ≈ 1` that wall-resolved LES
wants and the 30–100 of wall-*modeled* practice. Reaching `y+₁ ≈ 1` would need
`h_max/h_min ≈ 150–300`. Grading nearly doubles the points across `δ` (4.3 → 8.6),
which is the metric that matters here, but it does not change the regime.

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

All four configs therefore run at `dt = 1.0e-3` — for different reasons, some
CFL-limited and some margin-limited, not by copy-paste. If a Re=1e5 case proves
stable, raising it to 2.0e-3 halves its wall time and is the first lever to pull.

**Grading raises the CFL, and `dt` was deliberately left alone.** The smallest
element roughly halves, and the hot spot is the **step lip `(10, 0.75)`**, where
two walls meet (so wall distance → 0, giving `h_min`) *and* the gap jet runs at
`u ≈ 3`. Velocity and element size are helpfully correlated everywhere else —
both → 0 at walls, both large in the core — but not there. Predicted values, from
the `h_min` each mesh actually achieved:

| config | uniform CFL | predicted graded CFL |
|---|---|---|
| `bump2d.smoke` / `bump2d.Re500` | 0.192 | **~0.36** |
| `bump2d.Re1e5` | 0.096 | ~0.18 |
| `corner2d.smoke` / `corner2d.Re500` | 0.160 | **~0.32** |
| `corner2d.Re1e5` | 0.080 | ~0.16 |

Those are upper bounds — they assume peak velocity coincides with the smallest
element. **Run the smoke tests and read the reported CFL before launching
production.** If it lands near these numbers, halve `timestep_size`; the
Lax–Friedrichs dissipation may well carry ~0.36, but that is worth confirming for
an hour rather than assuming for a day.

For scale: our advective CFL is ~0.19 against OpenLB's 0.02, so at matched spatial
resolution the timestep is only about 3x larger. Landing in the same ballpark is
the sanity check that the mapping is right — it is not a speedup claim.

### Switching the inflow problem

Every config carries **both** `single_run` parameter blocks and selects between
them with `parameterized_problem/name`. Only the named block is read
(`ParameterizedProblem::SetSingleRun`), the same pattern
[`examples/stokes/array.8.yml`](../stokes/array.8.yml) uses.

**The inflow profiles are identical.** `channel_flow::ubdr` and
`backward_facing_step::ubdr` are the same Poiseuille parabola — verified to
machine precision (max difference 2.2e-16) — under different parameterizations:

| | `backward_facing_step` | `channel_flow` |
|---|---|---|
| formula | `u0 * 4/g^2 * (y-y0)(y1-y)`, `g = y1-y0` | `U * (1 - 4((y-x0)/L)^2)` |
| duct height | `y1 - y0` | `L` |
| duct centre | `(y0+y1)/2` | `x0` |
| peak velocity | `u0` | `U` |

**What differs is the initial condition.** `BackwardFacingStep` sets `ic_ptr`, so
the solver projects a Poiseuille field at `t = 0`. `ChannelFlow` sets none, and
`UnsteadyNSSolver::SetParameterizedProblem` then falls back to a zero IC — i.e.
**start from rest, as OpenLB does**. Which is preferable depends on the geometry:

- **bump2d**: `backward_facing_step` is the better fit and is the default. Its IC
  is divergence-free *and* consistent with the inflow BC at `t = 0` (the inlet
  spans the full channel height), so no ramp is needed. The only inconsistency is
  no-slip on the step faces, which relaxes within a few steps as a classic
  impulsive start.
- **corner2d**: that same IC is *not* consistent. With `y0=4, y1=5` it lays a jet
  across `x in [0,6]` at `y in [4,5]`, and that jet runs straight into step2's
  face at `x = 6`, violating no-penetration. It is resolved in the first implicit
  solve, but expect a large first-step pressure spike — watch the early CFL in the
  smoke test. `channel_flow` sidesteps it entirely and is closer to OpenLB.

**Why the mesh has to change too.** The two problems disagree about what boundary
tag 2 means, and a tag missing from a problem's `battr` leaves `bdr_type` unset and
trips `IsBdrTypeDefined()`, so no single tagging serves both:

| tag | `backward_facing_step` | `channel_flow` |
|---|---|---|
| 1 | inflow (Dirichlet) | wall (zero) |
| 2 | wall (zero) | **outflow (Neumann)** |
| 3 | outflow (Neumann) | wall (zero) |
| 4 | — | inflow (Dirichlet) |

So `generate_mesh.py --tags` emits two variants of every mesh — `.bfs` and
`.chan` — identical geometry, different physical group numbers. Switching means
changing `parameterized_problem/name` *and* `mesh/filename` together, which
`run_production.sh` does for you:

```bash
./run_production.sh corner2d.Re500 --problem channel_flow
```

One further asymmetry: `channel_flow` has no `amp`/`ky`/`freq` perturbation knobs,
so it cannot trip the flow at the inlet.

### What is *not* matched

**No turbulence model at Re = 1e5.** scaleupROM has no SGS model, so the Re=1e5
run is an under-resolved calculation stabilized by numerical dissipation: the
interior Lax–Friedrichs flux plays the role OpenLB gives to Smagorinsky. That is
an honest counterpart — the OpenLB run's molecular viscosity is negligible too
(`tau = 0.500018`) — but it is an implicit-LES comparison, so expect qualitative
agreement in separation and shedding behaviour, not pointwise agreement. In 2D
there is no vortex stretching and no real turbulent cascade on either side.

This is why all four configs set `full-discrete-galerkin: true`: the interior
Lax–Friedrichs integrator is only added when `full_dg` is on
(`SteadyNSSolver::BuildDomainOperators`). Without it the interior convection
operator is unstabilized central, which is fatal at a cell Reynolds number of
~3e3. At Re=500 the cell Reynolds number is only ~8, so those configs may safely
fall back to continuous Taylor–Hood (`full-discrete-galerkin: false`, which cuts
~760k DOF to ~310k for bump2d and ~700k to ~280k for corner2d) if the MUMPS
factorization does not fit in memory. **Do not do this at Re = 1e5** — it removes
the only stabilization in the scheme.

**No inflow ramp.** OpenLB starts from rest and ramps the inlet over the first 20
units via `PolynomialStartScale`. Neither problem here ramps: `channel_flow`
starts from rest but applies the inlet BC as a step in time, and
`backward_facing_step` starts from the fully-developed profile. Comparable
windows are therefore OpenLB `t in [20, 100]` against ours `t in [0, 80]`.

## Configuration notes

The problems are [`backward_facing_step` and `channel_flow`](../../src/parameterized_problem.cpp),
both used unmodified; see *Switching the inflow problem* above for the tag
conventions. The default, `backward_facing_step`, uses the same tag numbering as
[`utils/gmsh/bfs.geo`](../../utils/gmsh/bfs.geo): 1 = Dirichlet inflow, 2 = zero
(no-slip wall), 3 = Neumann outflow.

- `u0` (and `U`) is the Poiseuille **peak**, not the mean. `1.5` gives mean 1.0.
- `y0`/`y1` bound the inflow profile — `0.0`/`1.5` for bump2d, where the inlet
  spans the full channel, and `4.0`/`5.0` for corner2d's inlet duct.
- The viscosity that takes effect is `single_run/backward_facing_step/nu`, which
  overwrites `stokes/nu` during `SetParameterizedProblem`. Both are set to the
  same value in these configs to avoid confusion.
- `amp0`/`ky0`/`freq0`/`t_offset0` add a time-periodic transverse perturbation at
  the inlet. Off by default; try `amp0: 0.01, ky0: 1.0, freq0: 0.5` to trip the
  flow if it stays stubbornly steady.
- The Neumann outflow means `pres_dbc` is true, so the pressure-constant removal
  and complementary-flux paths stay off — matching OpenLB's pressure outflow BC.

## Meshes

Each geometry is declared in `generate_mesh.py` as a list of rectangular blocks.
The blocks no longer drive the mesher — the mesh is unstructured — but they still
define the fluid region, and the outline, area, solid holes and element-count
budget are all derived from them. Boundary curves are classified by the case's
inflow/outflow `x` locations, and the wall curves double as the size field's
`CurvesList`.

```bash
python3 generate_mesh.py --case bump2d   --h 0.0625  --check   # ~8,640 quads
python3 generate_mesh.py --case corner2d --h 0.03125 --check   # ~31,744 quads
python3 generate_mesh.py --case corner2d --tags chan --h 0.0625 --check
python3 generate_mesh.py --case bump2d   --ratio 1 --check      # uniform, for comparison
```

`--h` sets the DOF budget rather than the element size, and `--ratio` (default 3)
sets `h_max/h_min`. Since element count scales as `1/h_min²`, the generator tunes
`h_min` by secant iteration until the count is within 3% of the budget — normally
two passes. Achieved: `+1.3%` / `+0.7%` (bump2d) and `+2.4%` / `−1.0%` (corner2d).

`--check` re-reads the written `.msh` and asserts:

- **single element type** — 100% quads. This is load-bearing, not cosmetic:
  `SteadyNSSolver` requires one element type and builds its integration rules from
  a single `GeomType`, so one stray triangle is a real bug.
- **element count within 3%** of the budget — the fixed-DOF requirement.
- **boundary length per tag** equals the analytic perimeter split (bump2d
  54.5 = 1.5 inflow + 51.5 wall + 1.5 outflow; corner2d 32 = 1 + 30 + 1). This
  replaces the old edge-*count* invariant, which is meaningless once boundary
  spacing varies, and is strictly stronger.
- **grading direction** — wall-adjacent elements are at `h_min`, catching a
  transposed `SizeMin`/`SizeMax`, which is otherwise silent.
- **no interior refinement** — probes on the lines a structured mesher would have
  been forced to refine (bump2d `(2.5, 0.75)` and `(5, 1.15)`; corner2d
  `(3.5, 4.0)` and `(3.5, 1.0)`) must be at `h_max`, not `h_min`.
- conformity, no orphan nodes, exactly one surface tag (so the mesh loads as a
  single subdomain), inflow on `x_in` / outflow on `x_out` / walls on neither, and
  no element inside a solid region.

The script also prints `h_min`/`h_max`, the achieved count, the growth rate per
element, and the `timestep_size` implied by `h_min` alongside the uniform one — so
the CFL factor above can be re-derived at any resolution.

## If the smoke test disagrees

The load-bearing assumption in every `timestep_size` here is the `u_max` budget —
3.0 for bump2d, 2.5 for corner2d. The smoke tests report CFL every 20 steps. On
the **graded** meshes expect roughly **0.36** (bump2d) and **0.32** (corner2d);
the uniform-mesh baselines were 0.19 and 0.16. If the reading is near the graded
prediction, halve `timestep_size` in the production configs before launching. If
it comes in lower, the grading was nearly free and `dt` can stay.

For corner2d specifically, watch the **first few steps** rather than the settled
value: that is where `backward_facing_step`'s IC impulsively stops the inlet jet
at `x = 6`. A spike there is the signal to run that case with
`--problem channel_flow`, which starts from rest.

If a Re=1e5 run goes unstable, escalate in this order: halve `--dt`; switch to
`--problem channel_flow`; refine to `--h 0.03125`; add an inlet trip via
`amp0`/`ky0`/`freq0` (`backward_facing_step` only).
