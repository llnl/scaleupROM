# Incremental Parallel Solve Implementation Plan

## 1. Goal

Demonstrate that an existing scaleupROM problem can be assembled and solved
collectively on multiple MPI ranks, then establish a reusable distributed
iterative-solver path for Poisson, Stokes, Navier--Stokes, and other compatible
FOM systems.

The first milestone is not a topology redesign and is not tied to
advection--diffusion. It is a verified end-to-end parallel solve using the
infrastructure already present on the `parallel` branch.

Work in small, testable increments. Change production code only when a focused
run exposes a concrete failure or when the next end-to-end capability requires
the change. API cleanup, naming improvements, ownership abstractions, and storage
optimizations follow a working baseline rather than precede it.

This work also serves as a guided tour of the codebase. Every increment must
explain the relevant execution path, data ownership, configuration, and debugging
techniques well enough that a maintainer can reproduce and extend the work without
reverse-engineering the diff.

## 2. Success criteria

The initial capability is complete when one existing example:

1. starts under `mpirun` with at least two ranks;
2. has all ranks cooperate on the same problem;
3. assembles a distributed `HypreParMatrix`;
4. solves the system through an existing distributed solver path;
5. produces a solution and residual that agree with the one-rank reference; and
6. demonstrates that matrix-row ownership is distributed across ranks.

The preferred first example is parallel Poisson because it is the smallest
existing scalar path. Stokes is the second validation target. Navier--Stokes is
then enabled using the same shared solver and matrix infrastructure where its
operator permits it.

The preferred solver is a distributed iterative method. If that path is not
currently functional, distributed MUMPS is an acceptable temporary baseline.
The direct solve establishes assembly and ownership correctness; the iterative
path becomes the next milestone.

## 3. Working rules

### 3.1 Behavior first

Before editing production code:

1. state the exact command and expected behavior;
2. run the smallest relevant example or test;
3. record the observed failure or result;
4. make the smallest change that addresses that observation;
5. rerun the focused case; and
6. compare against the one-rank reference.

Do not introduce an abstraction solely because it may be useful later. A new API
must either be required by the next working solve increment or be deferred.

### 3.2 Small increments

Keep each production increment near 50 changed lines or fewer when practical.
Do not mix unrelated topology, matrix, solver, physics, I/O, and performance
changes. Tests and documentation may accompany the increment.

### 3.3 Evidence

For each run, record:

- exact command and input;
- MPI rank count;
- matrix global and local dimensions;
- solver and preconditioner;
- convergence status, iteration count, and final residual;
- comparison norm against the one-rank result; and
- any rank-local ownership diagnostics needed to prove distribution.

Do not describe a path as parallel merely because it launches under MPI. The
operator must be distributed and ranks must participate in the same solve.

### 3.4 Learning deliverable for every increment

Each increment must include a short maintainer walkthrough containing:

1. **How to run it:** build target, input file, launch command, rank count, and
   expected output files or terminal messages.
2. **Where execution goes:** the executable entry point and the important call
   chain through topology creation, mesh setup, assembly, matrix construction,
   solver setup, solve, and output.
3. **What each rank owns:** global versus rank-local objects, owned rows and
   subdomains, ghost/stored data, and which collectives exchange or reduce data.
4. **Why the code is structured this way:** the invariant or MFEM/Hypre
   requirement behind the implementation, not merely a description of lines
   changed.
5. **How to inspect it:** useful configuration options, logging, assertions,
   debugger breakpoints, and expected rank-local values.
6. **How to extend it:** the smallest relevant extension point for a new physics
   operator, linear solver, preconditioner, boundary condition, or input option.
7. **How it can fail:** common symptoms, likely causes, and the first diagnostic
   to run.

Keep these explanations close to the corresponding progress entry. When a topic
becomes stable and reusable, promote it into a durable maintainer document rather
than leaving it scattered across chronological notes.

### 3.5 Explain before changing

Before each production edit, state in plain language:

- the observed behavior;
- the relevant call path;
- the invariant that is violated or capability that is missing;
- why the proposed edit is the smallest appropriate change;
- what remains unchanged; and
- how the focused run distinguishes a real fix from an accidental pass.

Afterward, walk through the important diff and connect it back to the runtime
behavior. Refactoring must include the same explanation even when output is
intended to remain unchanged.

## 4. Stage 0: preserve and characterize the branch

### Increment 0A: inventory existing runnable paths

Inspect, without refactoring:

- `examples/poisson_parallel`;
- the Stokes example configured for Hypre assembly;
- `PoissonSolver::CreateHypreParMatrix()`;
- `StokesSolver::CreateHypreParMatrix()`;
- the serial and parallel MUMPS setup paths; and
- any existing Hypre Krylov configuration.

Identify the exact executable, input, and launch command that are intended to
exercise the current distributed path.

Also produce an initial code map showing:

```text
example main
  -> input/configuration parsing
  -> topology handler construction and load balancing
  -> physics solver construction
  -> domain and interface assembly
  -> HypreParMatrix construction
  -> linear solver setup and solve
  -> solution transfer, error calculation, and visualization
```

Use real function names and file links from the selected example. Mark which
steps execute independently on each rank and which invoke MPI collectives.

Completion: a short progress entry states what already exists, which example will
be the baseline, how to run it, and its actual call path.

### Increment 0B: establish one-rank reference

Run the selected example with one MPI rank without changing production code.
Record solution norm, residual, matrix dimensions, and output artifacts.

Explain how the input selects topology mode, matrix type, solver, boundary
conditions, and output. Identify the minimum parameters a maintainer normally
changes when defining a new problem.

Completion: the reference is reproducible and its configuration is understood.

### Increment 0C: run unchanged on two ranks

Run exactly the same problem with two ranks. This is the decisive first test of
the branch as it currently exists.

Classify the outcome:

- **Pass:** proceed to Stage 1 and add a regression test around the working path.
- **Assembly failure:** fix only the demonstrated distribution or indexing bug.
- **Solver failure:** retain the distributed matrix and use the existing
  distributed MUMPS path to separate matrix correctness from iterative-solver
  correctness.
- **Wrong result:** add the smallest matrix-action or solution-comparison test
  that localizes the disagreement before changing code.

Completion: the result is recorded with rank-local evidence, including failures.
The entry must explain what changes between a one-rank and two-rank run even when
the input is identical.

## 5. Stage 1: make the existing parallel solve reproducible

### Increment 1A: focused one-rank versus two-rank regression

Add the smallest practical automated or scripted comparison for the baseline
problem. Compare a physically meaningful solution norm and residual, not just
process exit status.

Completion: the test fails for any currently observed defect and passes after
the smallest necessary correction.

### Increment 1B: prove distributed ownership

Report each rank's owned row range and local matrix dimensions. Use existing
data and APIs first. Add diagnostics only if the current code cannot expose this
information.

Completion: each participating rank owns the expected non-overlapping row range,
and the ranges cover the global system.

Include one worked example mapping a global subdomain and global row to its owner,
rank-local index, stored neighbor data, and Hypre diagonal/off-diagonal block.

### Increment 1C: extend to four ranks

Run the same small problem on four ranks if the problem has enough subdomains.
Handle zero-owned-subdomain ranks only if this concrete run encounters them;
otherwise document the present restriction and defer generalization.

Completion: one-, two-, and four-rank results agree within a stated tolerance,
or a focused defect has been isolated for the next increment.

## 6. Stage 2: distributed iterative solver

Start this stage after distributed assembly has been validated independently,
using distributed MUMPS if necessary.

### Increment 2A: unpreconditioned Krylov correctness

On the Poisson baseline:

1. construct the Krylov solver with the topology communicator;
2. set the distributed `HypreParMatrix` as its operator;
3. use rank-local RHS and solution vectors matching the owned row interval;
4. set explicit tolerances and iteration limits; and
5. report convergence and the true final residual.

Use CG for an SPD Poisson operator. Use GMRES or FGMRES for nonsymmetric or
indefinite systems. Do not force one Krylov method onto every physics solver.

Completion: the unpreconditioned solve matches the validated direct solution on
one and two ranks. Performance is not yet a requirement.

The walkthrough must distinguish the matrix, solver, preconditioner, and vectors;
explain why their communicators and partitions must agree; and show where MFEM
wraps Hypre versus where scaleupROM supplies its own assembly logic.

### Increment 2B: first preconditioner

Add the smallest preconditioner appropriate to the baseline operator, preferably
BoomerAMG for Poisson. Record setup time, solve time, and iteration count.

Completion: the preconditioned solver converges reliably and retains numerical
agreement with the direct baseline.

Explain what the preconditioner approximates, why it is suitable for this
operator, and why the same choice may not work for Stokes or Navier--Stokes.

### Increment 2C: configuration and failure reporting

Expose only the options needed to select the working direct or iterative path,
tolerances, maximum iterations, and print level. A nonconverged solve must return
an explicit failure rather than silently continuing.

Completion: both validation-direct and iterative paths are selectable with a
small input change.

## 7. Stage 3: reuse across physics

### 3A: Stokes

Run the existing distributed Stokes assembly first with its known solver path.
Compare one- and two-rank solutions and residuals. Then select an iterative method
and preconditioner appropriate to its saddle-point block structure. Do not assume
the Poisson AMG configuration is sufficient.

Completion: distributed Stokes has a verified direct baseline and, if feasible,
a convergent iterative path.

Document the block unknowns and operator structure, how pressure null spaces or
constraints are handled, and how that structure changes solver selection compared
with scalar Poisson.

### 3B: Navier--Stokes

Treat each linearized Newton/Picard system as the distributed-solve unit. First
verify that the nonlinear algorithm calls a collective distributed linear solve
on every rank and that reductions use the solver communicator. Validate:

- the initial or Stokes-like linear solve;
- one nonlinear update;
- nonlinear residual reduction; and
- convergence against the one-rank result.

The nonsymmetric linearized operator will normally require GMRES/FGMRES and a
block or physics-aware preconditioner. This is a solver-specific extension, not
a reason to redesign topology first.

Completion: one small Navier--Stokes case completes collectively on two ranks and
matches the one-rank reference within stated nonlinear and linear tolerances.

Document the distinction between nonlinear iterations and their inner linear
solves, including which residual controls each loop and where convergence is
checked.

### 3C: other solvers

For another physics solver, reuse the same progression:

1. one-rank reference;
2. unchanged two-rank run;
3. distributed assembly validation;
4. distributed direct baseline if available;
5. operator-appropriate Krylov method;
6. operator-appropriate preconditioner; and
7. one-/two-/four-rank comparison.

Physics-specific assembly remains in the physics solver. Common solver setup may
be extracted only after at least two working users demonstrate the shared shape.

## 8. Deferred refactoring and optimization

The following are explicitly deferred until an end-to-end multi-rank solve works:

- renaming local counts to owned/stored/global terminology;
- adding a `DistributedLayout` abstraction;
- duplicating neighbor sets into new ghost metadata arrays;
- redesigning global/local index APIs;
- replacing complete ghost meshes with interface-only metadata;
- sparse containers for interface or ROM block arrays;
- eliminating every global metadata array;
- general load balancing beyond cases exercised by a working solver; and
- GPU/HIP enablement.

After the Poisson and Stokes paths work, characterize their repeated ownership,
mapping, communicator, and solver code. Refactor only duplication demonstrated by
those working implementations, preserving the regression tests established above.

## 9. ROM and workflow follow-up

Distributed ROM assembly and solve follow the same order:

1. validate one-rank ROM against its FOM/reference behavior;
2. run the existing Hypre ROM assembly unchanged on two ranks;
3. compare distributed matrix action with a gathered reference;
4. validate the existing distributed solve path;
5. introduce a distributed iterative solve if the reduced-system size warrants
   it; and
6. only then improve global reduced-vector storage or block containers.

Parallel sampling, snapshot I/O, communicator splitting, and concurrent training
jobs are separate milestones. They must not block proving that one physical
problem can be solved collectively by multiple ranks.

## 10. Immediate next actions

1. Do not build further work on the unstaged topology accessor/ghost-array
   refactor; retain or discard it separately from the solver milestone.
2. Identify the exact `examples/poisson_parallel` build target and smallest input.
3. provide the one-rank command and expected output;
4. provide the identical two-rank command;
5. inspect and record whether the run uses distributed MUMPS or a Hypre iterative
   solver; and
6. trace the example's real call path and explain its input options; and
7. make no production change until those runs reveal the actual first blocker.

## 11. Maintainer guide produced by this work

As knowledge stabilizes, consolidate the walkthroughs into a short guide with
these sections:

- building the relevant targets;
- running serial and MPI examples;
- reading and overriding input files;
- execution path from example to solver;
- topology, ownership, and ghost data;
- domain and interface assembly;
- MFEM sparse matrices versus Hypre distributed matrices;
- direct and iterative solver selection;
- interpreting residuals, errors, and output;
- adding a new problem or physics operator;
- adding a solver or preconditioner;
- testing one-rank versus multi-rank equivalence; and
- diagnosing common MPI failures.

This guide should describe verified behavior only. Open questions and temporary
limitations remain in the progress log until resolved.
