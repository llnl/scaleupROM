// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "hypre_utils.hpp"

#include <limits>
#include <vector>

namespace mfem
{

std::unique_ptr<HypreParMatrix> BuildHypreParMatrix(
   MPI_Comm comm, const SparseBlockMatrixLayout &layout,
   const SparseBlockMatrixDistribution &distribution,
   const SparseBlockMatrix &matrix)
{
   int rank = 0;
   int nranks = 0;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nranks);

   MFEM_VERIFY(distribution.NumRanks() == nranks,
               "Block distribution does not match the MPI communicator");
   MFEM_VERIFY(layout.NumBlockRows() == layout.NumBlockCols() &&
               layout.NumScalarRows() == layout.NumScalarCols(),
               "The Hypre solver adapter currently requires a square layout");
   MFEM_VERIFY(layout.RowOffsets().Size() == layout.ColOffsets().Size(),
               "Square Hypre layout must use matching block offsets");
   for (int i = 0; i < layout.RowOffsets().Size(); i++)
      MFEM_VERIFY(layout.RowOffsets()[i] == layout.ColOffsets()[i],
                  "Square Hypre layout must use matching block offsets");

   const int first_block = distribution.FirstBlock(rank);
   const int num_blocks = distribution.NumBlocks(rank);
   MFEM_VERIFY(matrix.NumStoredRows() == num_blocks,
               "Sparse matrix must contain every block row owned by this rank");
   for (int local_block_row = 0; local_block_row < num_blocks;
        local_block_row++)
      MFEM_VERIFY(matrix.GetGlobalRow(local_block_row) ==
                  first_block + local_block_row,
                  "Sparse matrix rows must follow rank-contiguous ownership order");

   size_t nnz = 0;
   for (int local_block_row = 0; local_block_row < num_blocks;
        local_block_row++)
      for (int entry = 0; entry < matrix.NumBlocksInRow(local_block_row); entry++)
         nnz += static_cast<size_t>(matrix.GetBlock(local_block_row, entry).
                                   NumNonZeroElems());
   MFEM_VERIFY(nnz <= static_cast<size_t>(std::numeric_limits<int>::max()),
               "Local nonzero count exceeds MFEM's integer CSR range");

   const int local_rows = distribution.NumScalars(rank);
   std::vector<int> row_starts(local_rows + 1);
   std::vector<HYPRE_BigInt> global_columns(nnz);
   std::vector<double> values(nnz);

   int local_scalar_row = 0;
   size_t output = 0;
   row_starts[0] = 0;
   for (int local_block_row = 0; local_block_row < num_blocks;
        local_block_row++)
   {
      const int global_block_row = matrix.GetGlobalRow(local_block_row);
      const int block_rows = layout.BlockRowSize(global_block_row);
      for (int block_row = 0; block_row < block_rows; block_row++)
      {
         for (int entry = 0; entry < matrix.NumBlocksInRow(local_block_row);
              entry++)
         {
            const int global_block_col =
               matrix.GetGlobalColumn(local_block_row, entry);
            const SparseMatrix &block = matrix.GetBlock(local_block_row, entry);
            const int *I = block.GetI();
            const int *J = block.GetJ();
            const double *data = block.GetData();
            const HYPRE_BigInt scalar_col_start =
               layout.BlockColStart(global_block_col);
            for (int j = I[block_row]; j < I[block_row + 1]; j++)
            {
               global_columns[output] = scalar_col_start + J[j];
               values[output] = data[j];
               output++;
            }
         }
         row_starts[++local_scalar_row] = static_cast<int>(output);
      }
   }
   MFEM_VERIFY(local_scalar_row == local_rows && output == nnz,
               "Sparse block-to-scalar CSR conversion produced inconsistent sizes");

   HYPRE_BigInt row_partition[2] = {
      distribution.FirstScalar(rank), distribution.EndScalar(rank)
   };
   HYPRE_BigInt col_partition[2] = {
      distribution.FirstScalar(rank), distribution.EndScalar(rank)
   };

   return std::unique_ptr<HypreParMatrix>(new HypreParMatrix(
      comm, local_rows, layout.NumScalarRows(), layout.NumScalarCols(),
      row_starts.data(), global_columns.data(), values.data(),
      row_partition, col_partition));
}

}
