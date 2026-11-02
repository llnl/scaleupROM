// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef LINALG_UTILS_HPP
#define LINALG_UTILS_HPP

#include "mfem.hpp"
#include "librom.h"
#include "etc.hpp"
#include <map>
#include <memory>
#include <vector>

using namespace mfem;

namespace CAROM
{

// Serial case interface between CAROM and MFEM
void ComputeCtAB(const Operator& A,
                 const CAROM::Matrix& B,  // Distributed matrix
                 const CAROM::Matrix& C,  // Distributed matrix
                 CAROM::Matrix& CtAB);    // Non-distributed (local) matrix

// Serial case set block of CAROM matrix.
void SetBlock(const CAROM::Matrix& blockMat,
               const int iStart, const int iEnd,
               const int jStart, const int jEnd,
               CAROM::Matrix& globalMat);

// Serial case set block of CAROM vector.
void SetBlock(const CAROM::Vector& blockVec,
               const int iStart, const int iEnd,
               CAROM::Vector& globalVec);

void CopyMatrix(const CAROM::Matrix &carom_mat,
                DenseMatrix &mfem_mat);

void PrintMatrix(const CAROM::Matrix &mat,
                 const std::string &filename);

// TODO: parallel version + hdf5 file format for larger basis dimension.
void PrintVector(const CAROM::Vector &vec,
                 const std::string &filename);

}

namespace mfem
{

struct MatrixBlocks
{
   int nrows;
   int ncols;
   Array2D<SparseMatrix *> blocks;

   MatrixBlocks() : nrows(0), ncols(0), blocks(0, 0) {}
   MatrixBlocks(const int nrows_, const int ncols_)
      : nrows(nrows_), ncols(ncols_), blocks(nrows_, ncols_)
   { blocks = NULL; }
   MatrixBlocks(const Array2D<SparseMatrix *> &mats)
   { *this = mats; }

   ~MatrixBlocks() { DeletePointers(blocks); }

   void SetSize(const int w, const int h)
   {
      for (int i = 0; i < nrows; i++)
         for (int j = 0; j < ncols; j++)
            if (blocks(i, j)) delete blocks(i, j);

      nrows = w; ncols = h;
      blocks.SetSize(w, h);
      blocks = NULL;
   }

   SparseMatrix*& operator()(const int i, const int j)
   {
      assert((i >= 0) && (i < nrows));
      assert((j >= 0) && (j < ncols));
      return blocks(i, j);
   }

   SparseMatrix const * operator()(const int i, const int j) const
   {
      assert((i >= 0) && (i < nrows));
      assert((j >= 0) && (j < ncols));
      return blocks(i, j);
   }

   // MatrixBlocks has the ownership. no copy assignment.
   MatrixBlocks& operator=(const Array2D<SparseMatrix *> &mats)
   {
printf("operator= const Array2D<SparseMatrix *>.\n");
      nrows = mats.NumRows();
      ncols = mats.NumCols();
      blocks = mats;

      return *this;
   }

   // Copy assignment.
   MatrixBlocks& operator=(const MatrixBlocks &matblock)
   {
printf("operator= const MatrixBlocks.\n");
      if (this == &matblock) return *this;
      DeletePointers(blocks);

      nrows = matblock.nrows;
      ncols = matblock.ncols;
      blocks.SetSize(nrows, ncols);
      blocks = NULL;

      for (int i = 0; i < nrows; i++)
         for (int j = 0; j < ncols; j++)
            if (matblock(i, j)) blocks(i, j) = new SparseMatrix(*matblock(i, j));

      return *this;
   }
};

/** Block-to-scalar numbering for a variable-size block matrix. */
class SparseBlockMatrixLayout
{
private:
   Array<int> global_row_offsets;
   Array<int> global_col_offsets;

   static void ValidateOffsets(const Array<int> &offsets,
                               const char *description);

public:
   explicit SparseBlockMatrixLayout(const Array<int> &global_offsets);
   SparseBlockMatrixLayout(const Array<int> &global_row_offsets_,
                           const Array<int> &global_col_offsets_);

   int NumBlockRows() const { return global_row_offsets.Size() - 1; }
   int NumBlockCols() const { return global_col_offsets.Size() - 1; }
   int NumScalarRows() const { return global_row_offsets.Last(); }
   int NumScalarCols() const { return global_col_offsets.Last(); }
   int BlockRowStart(const int global_block_row) const;
   int BlockColStart(const int global_block_col) const;
   int BlockRowSize(const int global_block_row) const;
   int BlockColSize(const int global_block_col) const;
   const Array<int> &RowOffsets() const { return global_row_offsets; }
   const Array<int> &ColOffsets() const { return global_col_offsets; }
};

class SparseBlockMatrixBuilder;

/** Immutable block-CSR data with globally indexed block rows and columns. */
class SparseBlockMatrix
{
private:
   Array<int> global_rows;
   Array<int> row_starts;
   Array<int> global_columns;
   std::vector<std::unique_ptr<SparseMatrix>> blocks;

   SparseBlockMatrix(const Array<int> &global_rows_,
                     const Array<int> &row_starts_,
                     const Array<int> &global_columns_,
                     std::vector<std::unique_ptr<SparseMatrix>> &&blocks_);
   friend class SparseBlockMatrixBuilder;

public:
   ~SparseBlockMatrix() = default;
   SparseBlockMatrix(const SparseBlockMatrix &) = delete;
   SparseBlockMatrix &operator=(const SparseBlockMatrix &) = delete;
   SparseBlockMatrix(SparseBlockMatrix &&) = default;
   SparseBlockMatrix &operator=(SparseBlockMatrix &&) = default;

   int NumStoredRows() const { return global_rows.Size(); }
   int GetGlobalRow(const int local_row) const;
   int NumBlocksInRow(const int local_row) const;
   int GetGlobalColumn(const int local_row, const int row_entry) const;
   const SparseMatrix &GetBlock(const int local_row, const int row_entry) const;
};

/** Mutable assembly object that produces an immutable SparseBlockMatrix. */
class SparseBlockMatrixBuilder
{
private:
   const SparseBlockMatrixLayout &layout;
   Array<int> global_rows;
   std::map<int, int> global_to_local_row;
   std::vector<std::map<int, std::unique_ptr<SparseMatrix>>> assembly_rows;
   bool consumed = false;

public:
   SparseBlockMatrixBuilder(const SparseBlockMatrixLayout &layout_,
                            const Array<int> &stored_global_rows);
   SparseBlockMatrixBuilder(const SparseBlockMatrixBuilder &) = delete;
   SparseBlockMatrixBuilder &operator=(const SparseBlockMatrixBuilder &) = delete;

   bool StoresGlobalRow(const int global_row) const;
   void AddBlock(const int global_row, const int global_col,
                 const SparseMatrix &block);
   SparseBlockMatrix Finalize();
};

void AddToBlockMatrix(const Array<int> &ridx, const Array<int> &cidx, const MatrixBlocks &mats, BlockMatrix &bmat);

void modifiedGramSchmidt(DenseMatrix& mat);

/* Orthonormalize mat over mat1 and itself. */
void Orthonormalize(DenseMatrix& mat1, DenseMatrix& mat);

// Compute Rt * A * P
void RtAP(DenseMatrix& R,
        const Operator& A,
        DenseMatrix& P,
        DenseMatrix& RAP);
DenseMatrix* DenseRtAP(DenseMatrix& R,
                       const Operator& A,
                       DenseMatrix& P);
SparseMatrix* SparseRtAP(DenseMatrix& R,
                        const Operator& A,
                        DenseMatrix& P);

void AddwRtAP(DenseMatrix& R,
            const Operator& A,
            DenseMatrix& P,
            DenseMatrix& RAP,
            const double w=1.0);
void AddwRtAP(DenseMatrix& R,
            const Operator& A,
            DenseMatrix& P,
            SparseMatrix& RAP,
            const double w=1.0);

template<typename T>
void PrintMatrix(const T &mat,
                const std::string &filename);

void PrintMatrix(const DenseMatrix &mat,
                 const std::string &filename);

// TODO: parallel version + hdf5 file format for larger basis dimension.
void PrintVector(const Vector &vec,
                 const std::string &filename);

SparseMatrix* ReadSparseMatrixFromHDF(const std::string filename);
void WriteSparseMatrixToHDF(const SparseMatrix* mat, const std::string filename);

// Matrix-Vector multiplication on the specified rows.
// Rows are not necessarily sorted and they can be duplicated.
void MultSubMatrix(const DenseMatrix &mat, const Array<int> &rows, const Vector &x, Vector &y);
void AddMultTransposeSubMatrix(const DenseMatrix &mat, const Array<int> &rows, const Vector &x, Vector &y);
void MultTransposeSubMatrix(const DenseMatrix &mat, const Array<int> &rows, const Vector &x, Vector &y);
// Currently needed only for DenseMatrix A.
// For a SparseMatrix A, we can utilize PartMult within this routine.
void AddSubMatrixRtAP(const DenseMatrix& R, const Array<int> &Rrows,
                      const DenseMatrix& A,
                      const DenseMatrix& P, const Array<int> &Prows,
                      DenseMatrix& RAP);
void AddSubMatrixRtAP(const DenseMatrix& R, const Array<int> &Rrows,
                      const DenseMatrix& A,
                      const DenseMatrix& P, const Array<int> &Prows,
                      SparseMatrix& RAP);

// DenseTensor is column major and i is the fastest index. 
// y_k = T_{ijk} * x_i * x_j
void TensorContract(const DenseTensor &tensor, const Vector &xi, const Vector &xj, Vector &yk);
// y_k += w * T_{ijk} * x_i * x_j
void TensorAddScaledContract(const DenseTensor &tensor, const double w, const Vector &xi, const Vector &xj, Vector &yk);
// Contracts along the last axis.
// M_{ij} = T_{ijk} * x_k
void TensorMult(const DenseTensor &tensor, const Vector &x, DenseMatrix &M);
// Contracts along the last axis.
// M_{ik} += T_{ijk} * x_j
void TensorAddMultOnJ(const DenseTensor &tensor, const Vector &x, DenseMatrix &M);
// Contracts along the axis (0 or 1) and add the multipled transpose.
// axis 0: M_{kj} += T_{ijk} * x_i
// axis 1: M_{ki} += T_{ijk} * x_j
void TensorAddMultTranspose(const DenseTensor &tensor, const Vector &x, const int axis, DenseMatrix &M);
// Contracts along the axis (0 or 1) and add the multipled transpose.
// axis 0: M_{kj} += w * T_{ijk} * x_i
// axis 1: M_{ki} += w * T_{ijk} * x_j
void TensorAddScaledMultTranspose(const DenseTensor &tensor, const double w, const Vector &x, const int axis, DenseMatrix &M);

class CGOptimizer : public OptimizationSolver
{
protected:
   const double golden_ratio = 0.5 * (1.0 + sqrt(5.0));
   const double Cr = 1.0 - 1. / golden_ratio;
   const double ib = 0.1;
   const double eps = 1.0e-14;

   const int N_mnbrak = 1e4;
   const int N_para = 1e4;

public:
   CGOptimizer() : OptimizationSolver() {}
   virtual ~CGOptimizer() {}

   virtual void SetOperator(const Operator &op) override
   {
      oper = &op;
      height = op.Height();
      width = op.Width();

      resvec.SetSize(height);
   }

   virtual void Mult(const Vector &rhs, Vector &sol) const;

private:
   mutable Vector resvec;
   mutable Vector sol1, g, xi, xi1, h;

   double Objective(const Vector &rhs, const Vector &sol) const;
   void Gradient(const Vector &rhs, const Vector &sol, Vector &grad) const;

   double Brent(const Vector &rhs, const Vector &xi, Vector &sol, double b0 = 1e-1, double lin_tol = 1e-2) const;
};

void GetBasisElement(const DenseMatrix &basis, const int col, const Array<int> vdofs,
                     Vector &basis_el, DofTransformation *dof_trans=NULL);

}

#endif
