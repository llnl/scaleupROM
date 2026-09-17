# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT

file(COPY bump2d.smoke.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)
file(COPY bump2d.Re500.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)
file(COPY bump2d.Re1e5.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)

file(COPY corner2d.smoke.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)
file(COPY corner2d.Re500.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)
file(COPY corner2d.Re1e5.yml DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)

file(COPY README.md DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)
file(COPY generate_mesh.py DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel)

file(COPY setup_highrechannel.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE)
file(COPY run_production.sh DESTINATION ${CMAKE_BINARY_DIR}/examples/highReChannel
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE)

# The meshes are not committed and not generated at build time: setup_highrechannel.sh
# produces them, so that the gmsh python module stays an example-time dependency
# rather than a build requirement.
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/examples/highReChannel/meshes)
