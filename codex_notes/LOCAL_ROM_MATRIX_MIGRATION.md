# Local ROM Matrix Migration

## Session handoff (2026-10-23)

This is the primary continuation document for the parallel component-ROM
matrix work. The work is intentionally being done in small increments with a
test at each boundary.

### Current state

- Increment 1, an ownership-agreement guard in
  `MultiBlockSolver::InitROMHandler`, is implemented.
- Increment 2, the generic `mfem::SparseBlockMatrix` container and its focused
  unit tests, is implemented.
- The production ROM assembly path has **not** been changed yet. It still
  assembles the replicated MFEM `BlockMatrix` and converts it to Hypre.
- There is no `LocalROMMatrix`, shadow assembly, production `Verify()` hook, or
  API taking two nullable matrix representations. Those prototype ideas were
  deliberately removed.
- No allocation-dependent parallel run is needed for Increment 2 because the
  new container is not wired into a solver yet.

The last validation completed successfully:

```bash
cd /usr/workspace/bergel1/scaleupROM/build
make -j 16 test_linalg_utils
cd test
./test_linalg_utils
```

All six tests passed, including the two `SparseBlockMatrixTest` cases, and
`git diff --check` passed. The build emits many pre-existing Intel compiler
warnings. The noninteractive test environment also emits harmless `stty` and
`basename` messages before GoogleTest starts.

### Why this work is relevant to the largest possible run

The target machine has roughly 112 CPU cores per node and approximately
256 GB of high-bandwidth memory per node. A failure or practical limit near 16
ranks means the application is not able to use the available machine.

The current global `BlockMatrix` contains an `N x N` dense table of block
pointers on every rank. The individual blocks are sparse and mostly absent,
but the pointer table is global and replicated. If the number of component
blocks per rank is fixed during weak scaling, then `N` grows with rank count
`P`, giving approximately:

```text
per-rank pointer metadata: O(N^2) = O(P^2)
aggregate pointer metadata: O(P * N^2) = O(P^3)
```

Thus, adding ranks also makes every rank's replicated allocation larger. This
can impose a problem-size ceiling even when the distributed numerical matrix
and solve would fit on the machine.

With sparse local block rows, each rank instead stores its owned rows and the
actual local/remote couplings in those rows. For a sparse component graph with
fixed work and coupling degree per rank, the intended weak-scaling behavior is
approximately constant matrix storage per rank and linear aggregate storage.
This change therefore directly supports larger component-ROM problems; it is
not merely a speed optimization.

It will not solve every scalability limitation. After the replicated matrix is
removed, likely next limits include global port traversal, other replicated
topology/operator metadata, basis/operator loading, and MUMPS factorization or
communication. Also, 112 MPI ranks per node may not be the eventual performance
optimum; fewer ranks plus threading may be faster. The code should nevertheless
be capable of using more ranks without replicated metadata exhausting memory.

### Exact next increment

Do **not** wire `SparseBlockMatrix` directly into Poisson assembly yet. First
introduce a shared ROM block partition/numbering object, with unit tests, that
answers:

- which global ROM block corresponds to `(global subdomain, variable)`;
- which rank owns a global ROM block;
- which global blocks are stored by the current rank;
- the global-to-local row mapping; and
- how domain-major and variable-major ordering behave.

Tests should cover one and multiple variables, uneven rank partitions,
domain-major ordering, variable-major ordering, and noncontiguous stored global
rows. This abstraction belongs in shared ROM infrastructure so Poisson, Stokes,
and Navier--Stokes use the same ownership rules.

A design decision must be made explicitly for variable-major ordering: owned
ROM rows may be noncontiguous in the existing global numbering, while a normal
Hypre row partition is contiguous. Possible solutions include changing the
distributed algebraic numbering to be rank-contiguous or maintaining an
explicit permutation. Do not silently assume Poisson's one-variable numbering.

After that object and its tests are stable, Increment 3 should add a Hypre
adapter consuming `SparseBlockMatrix` and compare its local diagonal/off-diagonal
CSR, column map, matrix action, and solution against the old path. Only after
equivalence tests pass should production parallel assembly stop allocating the
global `BlockMatrix`.

### Working-tree caution

The working tree contains unrelated untracked profiling scripts, build/editor
files, and core files. Preserve them. `src/unsteady_ns_solver.cpp` also appears
modified only because a missing final newline was added; it is unrelated to
this matrix work and should not be included accidentally when preparing a
focused commit.

## Goal

Replace the replicated global ROM `BlockMatrix` used by component-wise
parallel runs with rank-local block rows. Each rank should store and assemble
only:

- rows belonging to its owned component subdomains; and
- columns belonging either to an owned subdomain or to a remote subdomain
  coupled through a port.

The resulting local blocks will be used to construct Hypre's diagonal and
off-diagonal CSR matrices directly. The existing global matrix path will remain
available during development as a correctness reference and for serial/global
output where needed.

Poisson is the first validation problem because it has one variable and the
simplest block structure. The ownership model and local matrix container belong
in the shared `MultiBlockSolver`/`ROMHandler` infrastructure, not in
`PoissonSolver`. They must also support Stokes and Navier--Stokes, including
multiple variables, separate variable bases, nonlinear reduced operators, and
both domain-major and variable-major ROM ordering.

## Current execution path

The component ROM online path currently runs as follows:

```text
SingleRun
  -> MultiBlockSolver::AssembleROMMat
     -> allocate a global N x N MFEM BlockMatrix on every rank
     -> assemble local component rows and global port contributions
     -> MFEMROMHandler::CreateHypreParMatrix
        -> scan every global column for every owned row
        -> copy owned entries into Hypre diag/offd CSR matrices
  -> MFEMROMHandler::Solve
```

The final Hypre matrix and solve are distributed. The intermediate MFEM
`BlockMatrix` is not: its `Array2D` block-pointer table has `N x N` entries on
every rank even though the component graph is sparse.

The parallel Poisson FOM already follows the desired general pattern in
`PoissonSolver::CreateHypreParMatrix`: it holds owned and neighboring
subdomains, then fills Hypre diagonal and off-diagonal storage without first
constructing a replicated global block matrix.

## Invariants

For the first implementation target, component-wise Poisson with one basis per
subdomain:

1. ROM block `m` corresponds to global subdomain `m`.
2. Each global block row has exactly one owning MPI rank.
3. A rank may store remote block columns, but never remote block rows.
4. Both diagonal and off-diagonal interface contributions for an owned row are
   assembled by the row owner.
5. The union of owned row ranges covers every global block exactly once.

Topology and ROM ownership are currently computed independently by
`TopologyHandler::LoadBalance` and `MFEMROMHandler::LoadBalanceROMBlocks`.
Both use the same contiguous quotient/remainder partition for the initial
Poisson target, but this agreement was previously implicit.

## Increment 1: verify ownership agreement

Status: implemented.

`MultiBlockSolver::InitROMHandler` now checks, for component-wise ROMs without
separate variable bases, that:

- the ROM-owned block count equals the topology-owned subdomain count; and
- every topology-local subdomain maps to the corresponding global ROM block in
  the ROM handler's contiguous owned range.

This is deliberately behavior-preserving. It does not change assembly or
storage; it turns an assumption required by future local-row assembly into a
fail-fast invariant.

The check excludes separate-variable bases. Those have multiple ROM blocks per
subdomain and require an explicit ownership mapping that accounts for ROM block
ordering.

This exception is temporary and identifies the next shared design task; it is
not intended to make local assembly a Poisson-only feature.

### Validation

Build the project and run the existing component-wise workflow test. For a
parallel check, run the smallest available Poisson parallel configuration with
two ranks. Existing numerical results and matrix dimensions should remain
unchanged; any ownership disagreement should now fail during ROM handler
initialization rather than corrupting later assembly.

## Planned increments

## Increment 2: introduce generic sparse block-row storage

Status: core container and unit tests implemented; solver wiring intentionally
deferred.

The initial `LocalROMMatrix` shadow-assembly prototype was replaced with three
backend-neutral linear-algebra types. `SparseBlockMatrixLayout` owns the block
row/column offsets and translates block IDs to scalar starts and sizes.
`SparseBlockMatrixBuilder` accumulates contributions in sorted maps keyed by
global block-column ID. Its `Finalize()` method produces an immutable
`SparseBlockMatrix` containing only block-CSR-like row starts, global column
IDs, and contiguous row-oriented block storage. Production ROM assembly remains
on its original single-`BlockMatrix` interface.

The container supports rectangular matrices, variable block sizes,
noncontiguous stored global rows, and columns whose vector unknowns are owned by
another rank. It contains no MPI, ROM, topology, Hypre, or solver policy.
Allocation counts are checked before conversion to MFEM's local `int` index
type.

Focused GoogleTests in `test/test_linalg_utils.cpp` document and verify:

- variable-size block identity and dimensions;
- accumulation of repeated contributions to one block;
- deterministic row-oriented ordering after finalization;
- separation of rectangular layout metadata from finalized block data; and
- the row-ownership rule: a stored row may reference a global column whose
  input-vector unknowns are remote.

The focused test command is:

```bash
cd build/test
./test_linalg_utils --gtest_filter='SparseBlockMatrixTest.*'
```

All three focused tests and the complete `test_linalg_utils` executable pass. Since
this increment does not alter solver execution, no allocation-dependent
parallel run is required yet.

### Possible pattern-reuse follow-up

Keep the map-based `SparseBlockMatrixBuilder` for the initial production and
Hypre integration. Although the component and port topology should make the
block sparsity pattern predictable, that pattern still needs an initial
symbolic construction pass. The builder is the simpler starting point: it
discovers the block graph while accumulating contributions and produces the
block-CSR arrays during `Finalize()`.

If repeated matrix assembly or profiling later shows that rebuilding the block
graph is significant, introduce a separate immutable `SparseBlockMatrixPattern`
containing the stored global rows, row starts, and global block columns. Build
it once from topology (or capture it from the first builder result), then reuse
it to allocate subsequent matrices and to drive Hypre conversion without
global-column scans. Do not add this symbolic-pattern path before the dynamic
builder-based Hypre path is correct and measured.

### Ownership API follow-up

Before the local matrix becomes the Hypre input, shared ownership queries still
need to map `(global subdomain, variable)` to a global ROM block and owning rank
for both supported orderings. In particular, variable-major ordering can give a
rank noncontiguous global ROM block rows, while Hypre row partitions are
contiguous; that mismatch needs an explicit numbering/partition decision rather
than a Poisson-specific assumption.

## Increment 3: construct Hypre from local rows

Add a `CreateHypreParMatrix` overload that consumes local block rows. Compare
the resulting diagonal CSR, off-diagonal CSR, column map, matrix action, and
solution against the existing path on one and two ranks.

## Increment 4: stop allocating the global matrix in parallel

After equivalence tests pass, select local-row assembly for the distributed
Hypre path. Retain the old path for serial operation and explicit global matrix
output until those uses are handled separately.

## Increment 5: localize port traversal

Avoid scanning every global port on every rank. Traverse only ports touching an
owned row, using the topology handler's owned and neighbor information.

## Later solver extensions

- Stokes exercises multiple variables and an indefinite linear system.
- Steady Navier--Stokes reuses local linear block assembly and adds local
  nonlinear tensor or EQP operators plus a distributed nonlinear/Jacobian
  solve.
- Unsteady Navier--Stokes additionally reuses the same ownership for mass,
  advection, interface, and per-timestep linear systems.

These extensions should consume the shared ownership and local-row APIs. They
should not duplicate Poisson-specific Hypre conversion code.

## Validation checklist for every increment

- Build succeeds with assertions enabled.
- Existing serial component ROM workflow still passes.
- One-rank old and new matrix actions agree.
- Two-rank owned ranges are disjoint and cover the global rows.
- Two-rank old and new local Hypre matrices have matching dimensions and
  nonzero structure.
- The reduced solution agrees with the existing path within the established
  test tolerance.
- Profiling records show no unexpected regression before the old path is
  removed.
