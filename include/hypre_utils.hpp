// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HYPRE_UTILS_HPP
#define HYPRE_UTILS_HPP

#include "linalg_utils.hpp"

#include <memory>

namespace mfem
{

std::unique_ptr<HypreParMatrix> BuildHypreParMatrix(
   MPI_Comm comm, const SparseBlockMatrixLayout &layout,
   const SparseBlockMatrixDistribution &distribution,
   const SparseBlockMatrix &matrix);

}

#endif
