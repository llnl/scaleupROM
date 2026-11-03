// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include<gtest/gtest.h>
#include "linalg_utils.hpp"
#include "hypre_utils.hpp"
#include "etc.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
   SUCCEED();
}

TEST(linalg_test, MultSubMatrix)
{
   const int nrow = 15, ncol = 7, nsample = 4;
   DenseMatrix mat(nrow, ncol);
   Vector x(ncol), matx(nrow), y(nsample);
   Array<int> rows(nsample);

   for (int i = 0; i < nrow; i++)
      for (int j = 0; j < ncol; j++)
         mat(i, j) = UniformRandom();
   
   for (int j = 0; j < ncol; j++)
      x(j) = UniformRandom();

   for (int s = 0; s < nsample; s++)
      rows[s] = UniformRandom(0, nrow-1);

   y = 0.0;
   matx = 0.0;

   mat.Mult(x, matx);
   MultSubMatrix(mat, rows, x, y);

   for (int s = 0; s < nsample; s++)
      EXPECT_DOUBLE_EQ(y(s), matx(rows[s]));

   return;
}

TEST(linalg_test, AddMultTransposeSubMatrix)
{
   const int nrow = 15, ncol = 7, nsample = 4;
   DenseMatrix mat(nrow, ncol);
   Vector x(nsample), xtotal(nrow), matTx(ncol), y(ncol);
   Array<int> rows(nsample);

   for (int i = 0; i < nrow; i++)
      for (int j = 0; j < ncol; j++)
         mat(i, j) = UniformRandom();
   
   for (int j = 0; j < nsample; j++)
      x(j) = UniformRandom();

   for (int s = 0; s < nsample; s++)
      rows[s] = UniformRandom(0, nrow-1);

   xtotal = 0.0;
   for (int s = 0; s < nsample; s++)
      // there can be duplicated rows.
      xtotal(rows[s]) += x(s);

   y = 0.0;
   matTx = 0.0;

   mat.MultTranspose(xtotal, matTx);
   AddMultTransposeSubMatrix(mat, rows, x, y);

   for (int s = 0; s < ncol; s++)
      EXPECT_DOUBLE_EQ(y(s), matTx(s));

   return;
}

TEST(linalg_test, Orthonormalize)
{
   const double thre = 1.0e-14;
   const int nrow = 15, ncol1 = 7, ncol = 5;
   DenseMatrix mat1(nrow, ncol1);
   DenseMatrix mat(nrow, ncol);
   DenseMatrix test(ncol, ncol), test1(ncol1, ncol);

   for (int i = 0; i < nrow; i++)
   {
      for (int j = 0; j < ncol; j++)
         mat(i, j) = UniformRandom();

      for (int j = 0; j < ncol1; j++)
         mat1(i, j) = UniformRandom();
   }

   modifiedGramSchmidt(mat1);
   Vector scale(ncol1);
   for (int d = 0; d < ncol1; d++)
      scale(d) = UniformRandom();
   mat1.RightScaling(scale);

   Orthonormalize(mat1, mat);

   Vector tmp1, tmp2;
   for (int j = 0; j < ncol; j++)
   {
      mat.GetColumnReference(j, tmp1);
      test.GetColumnReference(j, tmp2);
      mat.MultTranspose(tmp1, tmp2);
   }
   for (int i = 0; i < ncol; i++)
      for (int j = 0; j < ncol; j++)
      {
         double val = (i == j) ? 1.0 : 0.0;
         EXPECT_NEAR(test(i, j), val, thre);
      };
   
   for (int j = 0; j < ncol; j++)
   {
      mat.GetColumnReference(j, tmp1);
      test1.GetColumnReference(j, tmp2);
      mat1.MultTranspose(tmp1, tmp2);
   }
   for (int i = 0; i < ncol1; i++)
      for (int j = 0; j < ncol; j++)
         EXPECT_NEAR(test1(i, j), 0.0, thre);

   return;
}

namespace
{

std::unique_ptr<SparseMatrix> MakeSparseMatrix(
   const int height, const int width,
   const std::initializer_list<std::tuple<int, int, double>> &entries)
{
   std::unique_ptr<SparseMatrix> matrix(new SparseMatrix(height, width));
   for (const auto &entry : entries)
      matrix->Set(std::get<0>(entry), std::get<1>(entry), std::get<2>(entry));
   return matrix;
}

}

TEST(SparseBlockMatrixTest, StoresVariableSizeBlocksByGlobalIdentity)
{
   // Three global blocks with scalar sizes 2, 1, and 2. This local matrix
   // stores global block rows 0 and 1 only.
   Array<int> offsets({0, 2, 3, 5});
   Array<int> stored_rows({0, 1});
   SparseBlockMatrixLayout layout(offsets);
   SparseBlockMatrixBuilder builder(layout, stored_rows);

   auto A00 = MakeSparseMatrix(2, 2, {{0, 0, 1.0}, {1, 1, 2.0}});
   auto A01 = MakeSparseMatrix(2, 1, {{0, 0, 3.0}, {1, 0, 4.0}});
   auto A11 = MakeSparseMatrix(1, 1, {{0, 0, 5.0}});

   builder.AddBlock(0, 0, *A00);
   builder.AddBlock(0, 1, *A01);
   builder.AddBlock(0, 1, *A01); // Contributions to one block accumulate.
   builder.AddBlock(1, 1, *A11);
   SparseBlockMatrix matrix = builder.Finalize();

   ASSERT_EQ(matrix.NumStoredRows(), 2);
   EXPECT_EQ(matrix.GetGlobalRow(0), 0);
   EXPECT_EQ(matrix.GetGlobalRow(1), 1);

   ASSERT_EQ(matrix.NumBlocksInRow(0), 2);
   EXPECT_EQ(matrix.GetGlobalColumn(0, 0), 0);
   EXPECT_EQ(matrix.GetGlobalColumn(0, 1), 1);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(0, 0)(1, 1), 2.0);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(0, 1)(0, 0), 6.0);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(0, 1)(1, 0), 8.0);

   ASSERT_EQ(matrix.NumBlocksInRow(1), 1);
   EXPECT_EQ(matrix.GetGlobalColumn(1, 0), 1);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(1, 0)(0, 0), 5.0);
}

TEST(SparseBlockMatrixTest, StoredRowsMayReferenceRemoteUnknownBlocks)
{
   // Think of rows 0 and 1 as rank 0's matrix rows. Global column block 2
   // represents input-vector unknowns owned by another rank. A12 still belongs
   // to rank 0 because matrix ownership is determined by its row.
   Array<int> offsets({0, 2, 4, 6, 8});
   Array<int> stored_rows({0, 1});
   SparseBlockMatrixLayout layout(offsets);
   SparseBlockMatrixBuilder builder(layout, stored_rows);
   auto A12 = MakeSparseMatrix(2, 2, {{0, 1, 7.0}, {1, 0, 8.0}});

   ASSERT_TRUE(builder.StoresGlobalRow(1));
   ASSERT_FALSE(builder.StoresGlobalRow(2));
   builder.AddBlock(1, 2, *A12);
   SparseBlockMatrix matrix = builder.Finalize();

   ASSERT_EQ(matrix.NumBlocksInRow(1), 1);
   EXPECT_EQ(matrix.GetGlobalColumn(1, 0), 2);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(1, 0)(0, 1), 7.0);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(1, 0)(1, 0), 8.0);
}

TEST(SparseBlockMatrixTest, KeepsLayoutSeparateFromFinalizedBlockData)
{
   Array<int> row_offsets({0, 2, 5});
   Array<int> col_offsets({0, 1, 5, 7});
   SparseBlockMatrixLayout layout(row_offsets, col_offsets);

   EXPECT_EQ(layout.NumBlockRows(), 2);
   EXPECT_EQ(layout.NumBlockCols(), 3);
   EXPECT_EQ(layout.NumScalarRows(), 5);
   EXPECT_EQ(layout.NumScalarCols(), 7);
   EXPECT_EQ(layout.BlockRowStart(1), 2);
   EXPECT_EQ(layout.BlockColStart(2), 5);
   EXPECT_EQ(layout.BlockRowSize(1), 3);
   EXPECT_EQ(layout.BlockColSize(2), 2);

   Array<int> stored_rows(1);
   stored_rows[0] = 1;
   SparseBlockMatrixBuilder builder(layout, stored_rows);
   auto A12 = MakeSparseMatrix(3, 2, {{2, 1, 9.0}});
   builder.AddBlock(1, 2, *A12);
   SparseBlockMatrix matrix = builder.Finalize();

   EXPECT_EQ(matrix.GetGlobalRow(0), 1);
   EXPECT_EQ(matrix.GetGlobalColumn(0, 0), 2);
   EXPECT_DOUBLE_EQ(matrix.GetBlock(0, 0)(2, 1), 9.0);
}

TEST(SparseBlockMatrixTest, DescribesRankContiguousBlockOwnership)
{
   Array<int> offsets({0, 2, 5, 6});
   SparseBlockMatrixLayout layout(offsets);
   Array<int> block_starts({0, 1, 3});
   SparseBlockMatrixDistribution distribution(layout, block_starts);

   EXPECT_EQ(distribution.NumRanks(), 2);
   EXPECT_EQ(distribution.FirstBlock(1), 1);
   EXPECT_EQ(distribution.EndBlock(1), 3);
   EXPECT_EQ(distribution.FirstScalarRow(1), 2);
   EXPECT_EQ(distribution.EndScalarRow(1), 6);
   EXPECT_EQ(distribution.NumScalarRows(1), 4);
   EXPECT_EQ(distribution.Owner(0), 0);
   EXPECT_EQ(distribution.Owner(2), 1);
   EXPECT_TRUE(distribution.OwnsBlock(1, 2));
   EXPECT_FALSE(distribution.OwnsBlock(0, 2));
}

TEST(SparseBlockMatrixTest, BuildsSerialHypreMatrixWithoutGlobalBlockScan)
{
   int nranks = 0;
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   if (nranks != 1) return;

   Array<int> offsets({0, 2, 3});
   SparseBlockMatrixLayout layout(offsets);
   Array<int> block_starts({0, 2});
   SparseBlockMatrixDistribution distribution(layout, block_starts);
   Array<int> stored_rows({0, 1});
   SparseBlockMatrixBuilder builder(layout, stored_rows);

   auto A00 = MakeSparseMatrix(2, 2, {{0, 0, 1.0}, {1, 1, 2.0}});
   auto A01 = MakeSparseMatrix(2, 1, {{0, 0, 3.0}, {1, 0, 4.0}});
   auto A10 = MakeSparseMatrix(1, 2, {{0, 0, 5.0}, {0, 1, 6.0}});
   auto A11 = MakeSparseMatrix(1, 1, {{0, 0, 7.0}});
   builder.AddBlock(0, 0, *A00);
   builder.AddBlock(0, 1, *A01);
   builder.AddBlock(1, 0, *A10);
   builder.AddBlock(1, 1, *A11);
   SparseBlockMatrix matrix = builder.Finalize();

   std::unique_ptr<HypreParMatrix> hypre =
      BuildHypreParMatrix(MPI_COMM_WORLD, layout, distribution, matrix);
   Vector x(3);
   x[0] = 1.0;
   x[1] = 2.0;
   x[2] = 3.0;
   Vector y(3);
   hypre->Mult(x, y);

   EXPECT_DOUBLE_EQ(y[0], 10.0);
   EXPECT_DOUBLE_EQ(y[1], 16.0);
   EXPECT_DOUBLE_EQ(y[2], 38.0);
}

TEST(SparseBlockMatrixTest, BuildsDistributedHypreMatrixFromLocalBlockRows)
{
   int rank = 0;
   int nranks = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   if (nranks != 2) return;

   // Interpret each scalar as a 1x1 block and distribute this global matrix by
   // block row across two ranks:
   //
   //                 [ 2  1 ] [ 1 ]   [ 4 ]
   //             A = [ 1  3 ] [ 2 ] = [ 7 ].
   //
   // Rank 0 stores block row [2 1], and rank 1 stores block row [1 3]. Each
   // rank therefore owns one scalar row and one entry of the distributed x and
   // y vectors, so its local Vector objects have size one. The off-diagonal
   // block on each rank references the other rank's x entry; HypreParMatrix::
   // Mult must communicate that remote value to produce local results 4 and 7.
   Array<int> offsets({0, 1, 2});
   SparseBlockMatrixLayout layout(offsets);
   Array<int> block_starts({0, 1, 2});
   SparseBlockMatrixDistribution distribution(layout, block_starts);
   Array<int> stored_rows(1);
   stored_rows[0] = rank;
   SparseBlockMatrixBuilder builder(layout, stored_rows);

   auto diagonal = MakeSparseMatrix(1, 1, {{0, 0, rank == 0 ? 2.0 : 3.0}});
   auto coupling = MakeSparseMatrix(1, 1, {{0, 0, 1.0}});
   builder.AddBlock(rank, rank, *diagonal);
   builder.AddBlock(rank, 1 - rank, *coupling);
   SparseBlockMatrix matrix = builder.Finalize();

   std::unique_ptr<HypreParMatrix> hypre =
      BuildHypreParMatrix(MPI_COMM_WORLD, layout, distribution, matrix);
   Vector x(1);
   x[0] = rank + 1.0;
   Vector y(1);
   hypre->Mult(x, y);

   EXPECT_DOUBLE_EQ(y[0], rank == 0 ? 4.0 : 7.0);
}

// TODO: add more tests from sketches/yaml_example.cpp.

int main(int argc, char* argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   MPI_Init(&argc, &argv);
   int result = RUN_ALL_TESTS();
   MPI_Finalize();
   return result;
}
