# Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: MIT

# Require CMake version 3.12 or later to expose modern features of
# CMake for maintainability. With additional work, it is possible to
# rewrite the build system generator for CMake 3.10+, 3.8+, 3.1+,
# 3.0+, or 2.8+, with increasing amounts of work as the minimum
# required version is decreased. Notably, CMake's find modules export
# targets (e.g., FindZLIB.cmake exports the ZLIB::ZLIB IMPORTED
# target) as of version 3.1, and CMake supports using the
# <PackageName>_ROOT variables to set the prefix search path for
# find_package(<PackageName>) as of version 3.12.
cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW) # Use <PackageName>_ROOT variables in find_package

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Even though this project is a pure C++ library, the C language must be enabled
# for MPI and HDF5 dependency checking. The Fortran language must be enabled
# for Fortran-C name mangling convention detection.
project(scaleupROM
  VERSION 1.0.0
  DESCRIPTION "Model reduction application to scaleUP flow simulations."
  HOMEPAGE_URL "https://lc.llnl.gov/gitlab/chung28/scaleupROM.git"
  LANGUAGES C CXX Fortran)

## Set a bunch of variables to generate a configure header
# Enable assertion checking if debug symbols generated
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEBUG_CHECK_ASSERTIONS "1")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

option(BUILD_SKETCHES "Build sketch files" OFF)

# # Define variables for use in generating a configure file
set(LIB_DIR "" CACHE PATH "User-defined path to the directory that contains all the prerequisites.")
set(MFEM_DIR "" CACHE PATH "User-defined path to mfem build.")
set(HYPRE_DIR "" CACHE PATH "User-defined path to hypre.")
set(METIS_DIR "" CACHE PATH "User-defined path to metis installation.")
set(PARMETIS_DIR "" CACHE PATH "User-defined path to parmetis installtion.")
set(MUMPS_DIR "" CACHE PATH "User-defined path to mumps.")
set(YAML_DIR "" CACHE PATH "User-defined path to yaml-cpp installation.")
set(LIBROM_DIR "" CACHE PATH "User-defined path to libROM source.")

# If LIB_DIR is specified and prerequisite path is not specified,
# then set a default path within LIB_DIR.
if (NOT (LIB_DIR STREQUAL ""))
  if (MFEM_DIR STREQUAL "")
    set(MFEM_DIR "${LIB_DIR}/mfem/build")
  endif()

  if (HYPRE_DIR STREQUAL "")
    set(HYPRE_DIR "${LIB_DIR}/hypre")
  endif()

  if (METIS_DIR STREQUAL "")
    set(METIS_DIR "${LIB_DIR}/metis-install")
  endif()

  if (PARMETIS_DIR STREQUAL "")
    set(PARMETIS_DIR "${LIB_DIR}/parmetis-install")
  endif()

  if (MUMPS_DIR STREQUAL "")
    set(MUMPS_DIR "${LIB_DIR}/mumps")
  endif()

  if (YAML_DIR STREQUAL "")
    set(YAML_DIR "${LIB_DIR}/yaml-cpp-install")
  endif()

  if (LIBROM_DIR STREQUAL "")
    set(LIBROM_DIR "${LIB_DIR}/libROM")
  endif()
endif()

# List minimum version requirements for dependencies where possible to make
# packaging easier later.
find_package(HDF5 1.8.0 REQUIRED)

find_package(BLAS 3.4.0 REQUIRED)
find_package(LAPACK 3.4.0 REQUIRED)

# In FindMPI.cmake, "CXX" refers to "the MPI C API [being] usable from C++"
find_package(MPI 1.2 REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)

find_package(Doxygen 1.8.5)

# Prefer the GoogleTest installation bundled under LIB_DIR.  On LC systems,
# other loaded modules (notably Python) can expose an ABI-incompatible GTest.
if (NOT (LIB_DIR STREQUAL "") AND NOT GTest_ROOT)
  set(GTest_ROOT "${LIB_DIR}")
endif()
find_package(GTest 1.6.0)

find_program(GMSH gmsh)

# MFEM is required.
find_library(MFEM mfem "${MFEM_DIR}" "${MFEM_DIR}/build" "${MFEM_DIR}/lib" "${MFEM_DIR}/build/lib" "$ENV{MFEM_DIR}" "$ENV{MFEM_DIR}/build" "$ENV{MFEM_DIR}/lib" "$ENV{MFEM_DIR}/build/lib")
find_library(HYPRE HYPRE "${HYPRE_DIR}/lib" "${HYPRE_DIR}/src/hypre/lib" "$ENV{HYPRE_DIR}/lib" "$ENV{HYPRE_DIR}/src/hypre/lib")
find_library(PARMETIS parmetis "${PARMETIS_DIR}/lib" "$ENV{PARMETIS_DIR}/lib")
find_library(METIS metis "${METIS_DIR}/lib" "$ENV{METIS_DIR}/lib")
find_library(MUMPS dmumps "${MUMPS_DIR}/lib" "${MUMPS_DIR}/build/local/lib" "$ENV{MUMPS_DIR}/lib" "$ENV{MUMPS_DIR}/build/local/lib")
find_path(MFEM_INCLUDES mfem.hpp "${MFEM_DIR}" "${MFEM_DIR}/build" "${MFEM_DIR}/include" "${MFEM_DIR}/build/include" "$ENV{MFEM_DIR}" "$ENV{MFEM_DIR}/build" "$ENV{MFEM_DIR}/include" "$ENV{MFEM_DIR}/build/include")
find_path(HYPRE_INCLUDES HYPRE.h "${HYPRE_DIR}/include" "${HYPRE_DIR}/src/hypre/include" "$ENV{HYPRE_DIR}/include" "$ENV{HYPRE_DIR}/src/hypre/include")
find_path(PARMETIS_INCLUDES metis.h "${METIS_DIR}/include" "${PARMETIS_DIR}/metis/include" "$ENV{METIS_DIR}/include"  "$ENV{PARMETIS_DIR}/metis/include")
find_path(MUMPS_INCLUDES dmumps_c.h "${MUMPS_DIR}/include" "${MUMPS_DIR}/build/local/include" "$ENV{MUMPS_DIR}/include" "$ENV{MUMPS_DIR}/build/local/include")

# yaml-cpp library
find_package(yaml-cpp REQUIRED HINTS "${YAML_DIR}" "$ENV{YAML_DIR}")
# find_library(YAML yaml-cpp HINTS "$ENV{YAML_DIR}/lib")
# find_path(YAML_INCLUDES yaml.h HINTS "$ENV{YAML_DIR}/include/yaml-cpp")

# libROM
find_library(LIBROM libROM.so HINTS "${LIBROM_DIR}/build/lib" "$ENV{LIBROM_DIR}/build/lib")
find_path(LIBROM_INCLUDES librom.h HINTS "${LIBROM_DIR}/lib" "$ENV{LIBROM_DIR}/lib")

include_directories(
  include
  ${MFEM_INCLUDES}
  ${HYPRE_INCLUDES}
  ${PARMETIS_INCLUDES}
  ${MUMPS_INCLUDES}
  ${HDF5_C_INCLUDE_DIRS}
  ${MPI_C_INCLUDE_DIRS}
  ${MFEM_C_INCLUDE_DIRS}
  ${LIBROM_INCLUDES}
)
link_libraries(
  ${MPI_C_LINK_FLAGS}
  ${MPI_C_LIBRARIES}
  MPI::MPI_C
  ${MPI_FORTRAN_LINK_FLAGS}
  ${MPI_FORTRAN_LIBRARIES}
  MPI::MPI_Fortran
  ${HDF5_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${MFEM}
  ${HYPRE}
  ${PARMETIS}
  ${METIS}
  ${MUMPS}
  yaml-cpp::yaml-cpp
  ${LIBROM}
)

set(scaleupROMObj_SOURCES

  include/etc.hpp
  src/etc.cpp

  include/hdf5_utils.hpp
  src/hdf5_utils.cpp

  include/block_smoother.hpp
  src/block_smoother.cpp

  include/rom_element_collection.hpp
  src/rom_element_collection.cpp

  include/unsteady_ns_solver.hpp
  src/unsteady_ns_solver.cpp

  include/topology_handler.hpp
  src/topology_handler.cpp

  include/component_topology_handler.hpp
  src/component_topology_handler.cpp

  include/interfaceinteg.hpp
  src/interfaceinteg.cpp

  include/dg_bilinear.hpp
  src/dg_bilinear.cpp

  include/dg_mixed_bilin.hpp
  src/dg_mixed_bilin.cpp

  include/dg_linear.hpp
  src/dg_linear.cpp

  include/nonlinear_integ.hpp
  src/nonlinear_integ.cpp

  include/hyperreduction_integ.hpp
  src/hyperreduction_integ.cpp

  include/rom_nonlinearform.hpp
  src/rom_nonlinearform.cpp

  include/interface_form.hpp
  src/interface_form.cpp

  include/rom_interfaceform.hpp
  src/rom_interfaceform.cpp

  include/multiblock_solver.hpp
  src/multiblock_solver.cpp

  include/poisson_solver.hpp
  src/poisson_solver.cpp

  include/linelast_solver.hpp
  src/linelast_solver.cpp

  include/stokes_solver.hpp
  src/stokes_solver.cpp

  include/steady_ns_solver.hpp
  src/steady_ns_solver.cpp
  src/steady_ns_solver_schwarz.cpp

  include/advdiff_solver.hpp
  src/advdiff_solver.cpp

  include/input_parser.hpp
  src/input_parser.cpp

  include/parameter.hpp
  src/parameter.cpp

  include/parameterized_problem.hpp
  src/parameterized_problem.cpp

  include/sample_generator.hpp
  src/sample_generator.cpp

  include/random_sample_generator.hpp
  src/random_sample_generator.cpp

  include/linalg_utils.hpp
  src/linalg_utils.cpp

  include/rom_handler.hpp
  src/rom_handler.cpp

  include/main_workflow.hpp
  src/main_workflow.cpp

)

add_library(scaleupROMObj OBJECT ${scaleupROMObj_SOURCES})
# target_link_libraries(scaleupROMObj
#   PUBLIC ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES} MPI::MPI_C ${MPI_FORTRAN_LINK_FLAGS} ${MPI_FORTRAN_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES}
#          ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MFEM} ${HYPRE} ${PARMETIS} ${METIS})
# target_include_directories(scaleupROMObj
#   PUBLIC ${MFEM_INCLUDES} ${HYPRE_INCLUDES} ${PARMETIS_INCLUDES} ${HDF5_C_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS} ${MFEM_C_INCLUDE_DIRS})

add_subdirectory(bin)
add_subdirectory(utils)
if (GTEST_FOUND)
  enable_testing()
  link_libraries(GTest::GTest)
  add_subdirectory(test)
endif()
add_subdirectory(examples)
if(BUILD_SKETCHES)
  add_subdirectory(sketches)
endif(BUILD_SKETCHES)

## Use the C++11 standard as an entire feature instead of
## enumerating individual compiler features for simplicity
#target_compile_features(ROM PRIVATE cxx_std_11)
