# Steady Navier--Stokes Parallel Status

## Current status

The steady Navier--Stokes FOM and ROM solve nonlinear algebraic systems using
Newton's method by default (L-BFGS is also available for the FOM). The current
implementation is suitable for serial use, but it is **not yet ready for a
correct or scalable distributed-memory solve**.

Launching the steady Navier--Stokes example with multiple MPI ranks does not by
itself establish spatial parallelism. Until the issues below are resolved, a
multi-rank run may perform redundant work, encounter incompatible matrix
layouts or communicators, or produce incorrectly coupled nonlinear systems.
Parallel timing results from this path should therefore not be treated as valid
strong- or weak-scaling data.

## Known gaps

- The FOM Newton Jacobian is constructed with `MPI_COMM_SELF`, and its direct
  MUMPS Jacobian solve also uses `MPI_COMM_SELF`. Ranks therefore solve their
  local/replicated systems independently rather than one distributed Newton
  system.
- Pressure mean removal uses a local vector sum and is marked `TODO` for
  parallelization. A distributed solve requires global reductions for pressure
  normalization.
- The ROM nonlinear solve creates its direct MUMPS solver on
  `MPI_COMM_WORLD`, while steady Navier--Stokes ROM Jacobians are currently
  wrapped with `MPI_COMM_SELF`. The nonlinear vector, Jacobian, and solver must
  use one consistent distributed layout and communicator.
- Tensor and empirical-quadrature (EQP) residual/Jacobian evaluation still need
  clear local-subdomain ownership and communication of interface contributions.
- Nonlinear EQP training and nonlinear artifact I/O still contain explicit
  parallelization TODOs.
- The distributed topology and ROM-matrix work validated for linear Poisson
  does not automatically distribute the steady Navier--Stokes nonlinear
  residual and Jacobian.

Relevant implementation locations include `src/steady_ns_solver.cpp`,
`include/steady_ns_solver.hpp`, and `src/rom_handler.cpp`.

## Recommended path to parallel support

1. Establish matching one-rank FOM and ROM reference solutions, including
   nonlinear residuals and Newton iteration histories.
2. Add a two-rank correctness test that compares the distributed solution and
   residual norm with the serial reference.
3. Partition nonlinear residual and Jacobian assembly by owned component
   blocks, with explicit exchange/accumulation of interface contributions.
4. Use collective reductions for pressure normalization and nonlinear norms.
5. Make Newton vectors, Jacobians, preconditioners, and linear solvers use the
   same `MPI_COMM_WORLD` ownership map.
6. Extend tensor/EQP loading, evaluation, and training to locally owned blocks
   and required neighbor data.
7. Only after correctness is established, measure end-to-end and nonlinear
   kernel strong/weak scaling separately.

## Minimum acceptance criteria

A steady Navier--Stokes path should be described as parallel only after it:

- agrees with a serial reference to a documented tolerance;
- gives rank-independent converged residuals and Newton iteration counts;
- assigns each component block to a well-defined owner rather than replicating
  the complete nonlinear calculation on every rank;
- uses compatible distributed communicators and matrix/vector layouts; and
- passes at least one multi-rank interface-coupling and pressure-normalization
  regression test.
