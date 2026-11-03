// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef ROM_HANDLER_HPP
#define ROM_HANDLER_HPP

#include "mfem.hpp"
#include "linalg/BasisGenerator.h"
#include "linalg/BasisReader.h"
#include "mfem/Utilities.hpp"
#include "topology_handler.hpp"
#include "linalg_utils.hpp"
#include "hdf5_utils.hpp"

namespace mfem
{

enum ROMBuildingLevel
{
   NONE,
   COMPONENT,
   GLOBAL,
   NUM_BLD_LVL
};

enum NonlinearHandling
{
   TENSOR,
   EQP,
   NUM_NLNHNDL
};

enum class ROMOrderBy
{
   VARIABLE,
   DOMAIN
};

const BasisTag GetBasisTagForComponent(const int &comp_idx, const TopologyHandler *topol_handler, const std::string var_name="");
const BasisTag GetBasisTag(const int &subdomain_index, const TopologyHandler *topol_handler, const std::string var_name="");

class ROMHandlerBase
{
protected:
// public:
   int numSub = -1;          // number of subdomains.
   int numSubLoc = -1;       // number of local subdomains.
   int num_var = -1;         // number of variables for which POD is performed.
   int num_rom_blocks = -1;  // number of ROM blocks for the global domain.
   int num_rom_blocks_local = -1;  // number of ROM blocks for the local domain.
   int num_rom_ref_blocks = -1;  // number of ROM reference component blocks.
   int num_rom_comp = -1;     // number of ROM reference components.
   std::vector<std::string> fom_var_names;          // dimension of each variable.
   Array<int> fom_var_offsets;
   Array<int> fom_num_vdofs;

   // rom options.
   bool save_sv = false;
   bool save_basis_visual = false;
   bool component_sampling = false;
   bool save_lspg_basis = false;
   ROMBuildingLevel save_operator = NUM_BLD_LVL;
   bool nonlinear_mode = false;
   bool separate_variable = false;
   NonlinearHandling nlin_handle = NUM_NLNHNDL;
   // ProjectionMode proj_mode = NUM_PROJMODE;

   // file names.
   std::string sample_dir;
   std::string sample_prefix;
   std::string basis_prefix;
   std::string operator_prefix;

   // topology handler
   TopologyHandler *topol_handler = NULL;

   // component rom variables.
   /*
      number of columns in the basis for a reference component.
      For i-th reference component and j-th variable,
         index = i * num_var + j
   */
   Array<int> num_ref_basis;
   /*
      number of rows in the basis for a reference component.
      For i-th reference component and j-th variable,
         index = i * num_var + j
   */
   Array<int> dim_ref_basis;
   Array<CAROM::Matrix*> carom_ref_basis;
   Array<int> rom_comp_block_offsets;
   std::vector<BasisTag> basis_tags;
   bool basis_loaded;
   bool operator_loaded;

   bool hypre_assemble = false;

   // MPI partitioning data
   Array<int> localSizes;
   Array<int> localNumBlocks;
   std::array<int,2> localBlocks;

   // domain rom variables.
   /*
      offset for the global domain ROM blocks.
      For i-th subdomain and j-th variable,
         ROMOrderBy::DOMAIN:   index = i * num_var + j
         ROMOrderBy::VARIABLE:     index = j * numSub + i
   */
   Array<int> num_basis;
   Array<int> rom_block_offsets;
   
   CAROM::Options* rom_options;
   CAROM::BasisGenerator *basis_generator;
   CAROM::BasisReader *basis_reader;

   int max_num_snapshots = 100;
   bool update_right_SV = false;
   bool incremental = false;

   ROMOrderBy ordering = ROMOrderBy::DOMAIN;

   mfem::BlockVector *reduced_rhs = NULL;
   mfem::BlockVector *reduced_sol = NULL;

   void ParseInputs();
public:
   ROMHandlerBase(TopologyHandler *input_topol, const Array<int> &input_var_offsets,
      const std::vector<std::string> &var_names, const bool separate_variable_basis);

   virtual ~ROMHandlerBase();

   // access
   const int GetNumSubdomains() { return numSub; }
   const int GetNumROMRefComps() { return num_rom_comp; }
   const int GetNumROMRefBlocks() { return num_rom_ref_blocks; }
   const int GetRefNumBasis(const int &basis_idx) { return num_ref_basis[basis_idx]; }
   const ROMBuildingLevel GetBuildingLevel() { return save_operator; }
   const bool BasisLoaded() { return basis_loaded; }
   const bool OperatorLoaded() { return operator_loaded; }
   const bool SeparateVariable() { return separate_variable; }
   const std::string GetOperatorPrefix() { return operator_prefix; }
   const std::string GetBasisPrefix() { return basis_prefix; }
   const BasisTag GetRefBasisTag(const int ref_idx) { return basis_tags[ref_idx]; }
   const Array<int>* GetBlockOffsets() { return &rom_block_offsets; }
   virtual SparseMatrix* GetOperator() = 0;
   const bool GetNonlinearMode() { return nonlinear_mode; }
   void SetNonlinearMode(const bool nl_mode)
   {
      if (nlin_handle == NonlinearHandling::NUM_NLNHNDL)
         mfem_error("ROMHandler::SetNonlinearMode - nonlinear handling is not set!\n");
      nonlinear_mode = nl_mode;
   }
   const NonlinearHandling GetNonlinearHandling() { return nlin_handle; }
   const ROMOrderBy GetOrdering() { return ordering; }
   const int GetBlockIndex(const int m, const int v=-1);
   void GetDomainAndVariableIndex(const int &rom_block_index, int &m, int &v);

   bool HypreAssemble() const { return hypre_assemble; }

   std::array<int,2> GetLocalBlocks() { return localBlocks; }
   void GetLocalNumBlocks(Array<int> &lnb) { lnb = localNumBlocks; }

   mfem::BlockVector* GetReducedSolution() { return reduced_sol; }
   mfem::BlockVector* GetReducedRHS() { return reduced_rhs; }

   /* parse inputs for supremizer. only for Stokes/SteadyNS Solver. */
   void ParseSupremizerInput(Array<int> &num_ref_supreme, Array<int> &num_supreme);

   virtual void LoadReducedBasis();

   int GetRefIndexForSubdomain(const int &subdomain_index);
   virtual void GetReferenceBasis(const int &basis_index, DenseMatrix* &basis) = 0;
   virtual void GetDomainBasis(const int &basis_index, DenseMatrix* &basis) = 0;
   virtual void SetBlockSizes();

   // P_i^T * mat * P_j
   virtual SparseMatrix* ProjectToRefBasis(const int &i, const int &j, const Operator *mat) = 0;
   virtual SparseMatrix* ProjectToDomainBasis(const int &i, const int &j, const Operator *mat) = 0;
   virtual void ProjectToRefBasis(const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectToDomainBasis(const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;

   virtual void ProjectComponentToRefBasis(const int &c, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectComponentToDomainBasis(const int &m, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectInterfaceToRefBasis(const int &c1, const int &c2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectInterfaceToDomainBasis(const int &m1, const int &m2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectVariableToDomainBasis(const int &vi, const int &vj, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectGlobalToDomainBasis(const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats) = 0;
   virtual void ProjectOperatorOnReducedBasis(const Array2D<Operator*> &mats) = 0;

   virtual void ProjectToRefBasis(const int &i, const Vector &vec, Vector &rom_vec) = 0;
   virtual void ProjectToDomainBasis(const int &i, const Vector &vec, Vector &rom_vec) = 0;
   virtual void ProjectGlobalToDomainBasis(const BlockVector* vec, mfem::BlockVector*& rom_vec) = 0;
   virtual void ProjectRHSOnReducedBasis(const BlockVector* RHS) = 0;

   virtual void LiftUpFromRefBasis(const int &i, const Vector &rom_vec, Vector &vec) = 0;
   virtual void LiftUpFromDomainBasis(const int &i, const Vector &rom_vec, Vector &vec) = 0;
   virtual void LiftUpGlobal(const BlockVector &rom_vec, BlockVector &vec) = 0;

   virtual void Solve(Vector &rhs, Vector &sol) = 0;
   virtual void Solve(BlockVector* U) = 0;
   virtual void NonlinearSolve(Operator &oper, BlockVector* U, Solver *prec=NULL) = 0;   

   virtual void SaveOperator(const std::string filename) = 0;
   virtual void LoadOperatorFromFile(const std::string filename) = 0;
   virtual void SetRomMat(BlockMatrix *input_mat, const bool init_direct_solver=true) = 0;
   virtual void SaveRomSystem(const std::string &input_prefix, const std::string type="mm") = 0;

   virtual void SaveBasisVisualization(const Array<FiniteElementSpace *> &fes, const std::vector<std::string> &var_names) = 0;

   virtual void SaveReducedSolution(const std::string &filename) = 0;
   virtual void SaveReducedRHS(const std::string &filename) = 0;

   virtual void AppendReferenceBasis(const int &idx, const DenseMatrix &mat) = 0;

   virtual void CreateHypreParMatrix(BlockMatrix *input_mat, int rank, int nproc) = 0;
   virtual void CreateHypreParMatrix(
      const SparseBlockMatrixLayout &layout,
      const SparseBlockMatrixDistribution &distribution,
      const SparseBlockMatrix &matrix) = 0;
   virtual void LoadBalanceROMBlocks(int rank, int nproc) = 0;
};

class MFEMROMHandler : public ROMHandlerBase
{
protected:
   // type for linear system solver.
   enum SolverType
   {
      DIRECT,
      CG,
      MINRES,
      GMRES,
      NUM_SOLVERTYPE
   } linsol_type;
   MUMPSSolver::MatType mat_type;

   // component rom variables.
   Array<DenseMatrix*> ref_basis;

   // domain rom variables.
   Array<DenseMatrix*> dom_basis; // This is only the pointers to ref_basis. no need of deleting.

   BlockMatrix *romMat = NULL;
   SparseMatrix *romMat_mono = NULL;

   // hypre matrix data
   SparseMatrix *hdiag = NULL;
   SparseMatrix *hoffd = NULL;
   HYPRE_BigInt *cmap = NULL;

   // variables needed for direct solve
   HYPRE_BigInt sys_glob_size;
   HYPRE_BigInt sys_row_starts[2];
   HypreParMatrix *romMat_hypre = NULL;
   MUMPSSolver *mumps = NULL;

   Vector reduced_rhs_hypre;
   HYPRE_BigInt hypre_start;

public:
   MFEMROMHandler(TopologyHandler *input_topol, const Array<int> &input_var_offsets,
      const std::vector<std::string> &var_names, const bool separate_variable_basis);

   virtual ~MFEMROMHandler();

   virtual SparseMatrix* GetOperator() override
   { assert(romMat_mono); return romMat_mono; }
   
   virtual void LoadReducedBasis();
   virtual void GetReferenceBasis(const int &basis_index, DenseMatrix* &basis) override;
   virtual void GetDomainBasis(const int &basis_index, DenseMatrix* &basis);

   // P_i^T * mat * P_j
   virtual SparseMatrix* ProjectToRefBasis(const int &i, const int &j, const Operator *mat);
   virtual SparseMatrix* ProjectToDomainBasis(const int &i, const int &j, const Operator *mat);
   virtual void ProjectToRefBasis(const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectToDomainBasis(const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);

   virtual void ProjectComponentToRefBasis(const int &c, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectComponentToDomainBasis(const int &m, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectInterfaceToRefBasis(const int &c1, const int &c2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectInterfaceToDomainBasis(const int &m1, const int &m2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectVariableToDomainBasis(const int &vi, const int &vj, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectGlobalToDomainBasis(const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats);
   virtual void ProjectOperatorOnReducedBasis(const Array2D<Operator*> &mats);

   virtual void ProjectToRefBasis(const int &i, const Vector &vec, Vector &rom_vec);
   virtual void ProjectToDomainBasis(const int &i, const Vector &vec, Vector &rom_vec);
   virtual void ProjectGlobalToDomainBasis(const BlockVector* vec, mfem::BlockVector*& rom_vec);
   virtual void ProjectRHSOnReducedBasis(const BlockVector* RHS) override
   { ProjectGlobalToDomainBasis(RHS, reduced_rhs); }

   virtual void LiftUpFromRefBasis(const int &i, const Vector &rom_vec, Vector &vec);
   virtual void LiftUpFromDomainBasis(const int &i, const Vector &rom_vec, Vector &vec);
   virtual void LiftUpGlobal(const BlockVector &rom_vec, BlockVector &vec);
   
   void Solve(Vector &rhs, Vector &sol) override;
   void Solve(BlockVector* U) override;
   void NonlinearSolve(Operator &oper, BlockVector* U, Solver *prec=NULL) override;

   virtual void SaveOperator(const std::string input_prefix="");
   virtual void LoadOperatorFromFile(const std::string input_prefix="");
   virtual void SetRomMat(BlockMatrix *input_mat, const bool init_direct_solver=true);
   virtual void SaveRomSystem(const std::string &input_prefix, const std::string type="mm");

   virtual void SaveBasisVisualization(const Array<FiniteElementSpace *> &fes, const std::vector<std::string> &var_names);

   virtual void SaveReducedSolution(const std::string &filename) override
   { PrintVector(*reduced_sol, filename); }
   virtual void SaveReducedRHS(const std::string &filename) override
   { PrintVector(*reduced_rhs, filename); }

   virtual void AppendReferenceBasis(const int &idx, const DenseMatrix &mat);

   void CreateHypreParMatrix(BlockMatrix *input_mat, int rank, int nproc) override;
   void CreateHypreParMatrix(
      const SparseBlockMatrixLayout &layout,
      const SparseBlockMatrixDistribution &distribution,
      const SparseBlockMatrix &matrix) override;
   void LoadBalanceROMBlocks(int rank, int nproc) override;

private:
   IterativeSolver* SetIterativeSolver(const MFEMROMHandler::SolverType &linsol_type_, const std::string &prec_type);
   // void GetBlockSparsity(const SparseMatrix *mat, const Array<int> &block_offsets, Array2D<bool> &mat_zero_blocks);
   // bool CheckZeroBlock(const DenseMatrix &mat);
   void SetupDirectSolver();
};


} // namespace mfem

#endif
