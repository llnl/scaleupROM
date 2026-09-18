# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT

file(COPY square-circle.geo DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY square.mesh DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY square.tri.mesh DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY test.multi_comp.h5 DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY box-channel.2x2.periodic.h5 DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY box-channel.1x2.h5 DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY test.component.yml DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY stokes.component.yml DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY steadyns.lf.yml DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY steadyns.interface_eqp.yml DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)
file(COPY usns.periodic.yml DESTINATION ${CMAKE_BINARY_DIR}/test/gmsh/)

ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_BINARY_DIR}/test/gmsh/square-circle.msh
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/gmsh"
    COMMAND "${GMSH}" square-circle.geo -2
    COMMENT "Generating test/gmsh/square-circle.msh ..."
)

ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_BINARY_DIR}/test/gmsh/square-circle.msh.mfem
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/test/gmsh"
    COMMAND ${CMAKE_BINARY_DIR}/utils/gmsh2mfem -m square-circle.msh
    DEPENDS ${CMAKE_BINARY_DIR}/utils/gmsh2mfem ${CMAKE_BINARY_DIR}/test/gmsh/square-circle.msh
    COMMENT "Generating test/gmsh/square-circle.msh.mfem ..."
)

add_executable(multi_comp_dd_mms multi_comp_dd_mms.cpp
                $<TARGET_OBJECTS:scaleupROMObj>
                $<TARGET_OBJECTS:mmsSuiteObj>
                ${CMAKE_BINARY_DIR}/test/gmsh/square-circle.msh.mfem
                )

add_executable(stokes_multi_comp_dd_mms stokes_multi_comp_dd_mms.cpp
                $<TARGET_OBJECTS:scaleupROMObj>
                $<TARGET_OBJECTS:mmsSuiteObj>
                )

add_executable(test_multi_comp_workflow test_multi_comp_workflow.cpp
                $<TARGET_OBJECTS:scaleupROMObj>
                )
