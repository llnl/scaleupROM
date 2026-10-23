# Copyright 2024 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT

file(COPY poisson.sampling.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY array.8.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY benchmark.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)

file(COPY generate_configs.py DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY setup_poisson.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY benchmark.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY clear_results.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY run_workflow.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY generate_benchmark.py DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY generate_structured_benchmark.py DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY plot_benchmarks.py DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY profile_benchmarks.py DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)
file(COPY benchmarks DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel)

file(COPY ../stokes/meshes/square.o3.mesh DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel/meshes)
file(COPY ../stokes/meshes/square-circle.msh.mfem DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel/meshes)
file(COPY ../stokes/meshes/square-triangle.msh.mfem DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel/meshes)
file(COPY ../stokes/meshes/square-star.msh.mfem DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel/meshes)
file(COPY ../stokes/meshes/square-square.msh.mfem DESTINATION ${CMAKE_BINARY_DIR}/examples/poisson_parallel/meshes)
