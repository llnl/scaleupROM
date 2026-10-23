// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "etc.hpp"
#include "input_parser.hpp"
#include "rom_handler.hpp"
#include "hdf5_utils.hpp"
#include "block_smoother.hpp"
#include "utils/mpi_utils.h"  // this is from libROM/utils.
// #include <cmath>
// #include <algorithm>

using namespace std;

namespace mfem
{

const BasisTag GetBasisTagForComponent(
   const int &comp_idx, const TopologyHandler *topol_handler, const std::string var_name)
{
   BasisTag tag;
   tag.comp = topol_handler->GetComponentName(comp_idx);
   if (var_name != "")
      tag.var = var_name;
   return tag;
}

const BasisTag GetBasisTag(
   const int &subdomain_index, const TopologyHandler *topol_handler, const std::string var_name)
{
   int c_type = topol_handler->GetMeshType(subdomain_index);
   return GetBasisTagForComponent(c_type, topol_handler, var_name);
}

ROMHandlerBase::ROMHandlerBase(
   TopologyHandler *input_topol, const Array<int> &input_var_offsets,
   const std::vector<std::string> &var_names, const bool separate_variable_basis)
   : topol_handler(input_topol),
     numSub(input_topol->GetNumSubdomains()),
     numSubLoc(input_topol->GetNumLocalSubdomains()),
     fom_var_names(var_names),
     fom_var_offsets(input_var_offsets),
     basis_tags(0),
     basis_loaded(false),
     operator_loaded(false),
     separate_variable(separate_variable_basis)
{
   num_var = fom_var_names.size();

   fom_num_vdofs.SetSize(numSubLoc);
   for (int m = 0, idx = 0; m < numSubLoc; m++)
   {
      fom_num_vdofs[m] = 0;
      for (int v = 0; v < num_var; v++, idx++)
         fom_num_vdofs[m] += fom_var_offsets[idx+1] - fom_var_offsets[idx];
   }

   ParseInputs();

   SetBlockSizes();
}

ROMHandlerBase::~ROMHandlerBase()
{
   DeletePointers(carom_ref_basis);
   delete reduced_rhs;
   delete reduced_sol;
}

void ROMHandlerBase::ParseInputs()
{
   assert(numSub > 0);
   assert(topol_handler != NULL);

   basis_prefix = config.GetOption<std::string>("basis/prefix", "basis");

   std::string save_op_str = config.GetOption<std::string>("model_reduction/save_operator/level", "none");
   if (save_op_str == "none")             save_operator = ROMBuildingLevel::NONE;
   else if (save_op_str == "global")      save_operator = ROMBuildingLevel::GLOBAL;
   else if (save_op_str == "component")   save_operator = ROMBuildingLevel::COMPONENT;
   else
      mfem_error("Unknown ROM building level!\n");

   if (save_operator != ROMBuildingLevel::NONE)
      operator_prefix = config.GetRequiredOption<std::string>("model_reduction/save_operator/prefix");

   num_rom_blocks = numSub;
   num_rom_blocks_local = num_rom_blocks;
   num_rom_ref_blocks = topol_handler->GetNumComponents();

   num_rom_comp = num_rom_ref_blocks;

   if (separate_variable)
   {
      num_rom_ref_blocks *= num_var;
      num_rom_blocks *= num_var;
   }
   
   const int num_ref_basis_default = config.GetOption<int>("basis/number_of_basis", -1);
   num_ref_basis.SetSize(num_rom_ref_blocks);
   num_ref_basis = num_ref_basis_default;
   assert(num_ref_basis.Size() > 0);
   dim_ref_basis.SetSize(num_rom_ref_blocks);

   YAML::Node basis_list = config.FindNode("basis/tags");
   if ((!basis_list) && (num_ref_basis_default <= 0))
      mfem_error("ROMHandler - cannot find the basis tag list, nor default number of basis is not set!\n");

   for (int b = 0; b < num_rom_ref_blocks; b++)
   {
      /* determine basis tag */
      BasisTag basis_tag;
      if (separate_variable)
         basis_tag = GetBasisTagForComponent(b / num_var, topol_handler, fom_var_names[b % num_var]);
      else
         basis_tag = GetBasisTagForComponent(b, topol_handler);

      /* determine number of basis */
      num_ref_basis[b] = config.LookUpFromDict("name", basis_tag.print(), "number_of_basis", num_ref_basis_default, basis_list);

      if (num_ref_basis[b] < 0)
      {
         printf("Cannot find the number of basis for %s!\n", basis_tag.print().c_str());
         mfem_error("Or specify the default number of basis!\n");
      }

      /* parse the dimension of basis */
      int midx = -1;
      int vidx = (separate_variable) ? b % num_var : 0;
      for (int m = 0; m < numSubLoc; m++)
         if (topol_handler->GetMeshType(m) == (separate_variable ? b / num_var : b))
         {
            midx = m;
            break;
         }
      assert(midx >= 0);

      const int idx = (separate_variable) ? midx * num_var + vidx : midx;
      if (separate_variable)
         dim_ref_basis[b] = fom_var_offsets[idx + 1] - fom_var_offsets[idx];
      else
         dim_ref_basis[b] = fom_num_vdofs[idx];
   }  // for (int b = 0; b < num_rom_ref_blocks; b++)
   
   for (int k = 0; k < num_ref_basis.Size(); k++)
   {
      assert(num_ref_basis[k] > 0);
      assert(dim_ref_basis[k] > 0);
   }

   max_num_snapshots = config.GetOption<int>("sample_generation/maximum_number_of_snapshots", 100);
   update_right_SV = config.GetOption<bool>("basis/svd/update_right_sv", false);

   save_sv = config.GetOption<bool>("basis/svd/save_spectrum", false);

   save_basis_visual = config.GetOption<bool>("basis/visualization/enabled", false);

   std::string nlin_handle_str = config.GetOption<std::string>("model_reduction/nonlinear_handling", "none");
   if (nlin_handle_str == "tensor")    nlin_handle = NonlinearHandling::TENSOR;
   else if (nlin_handle_str == "eqp")  nlin_handle = NonlinearHandling::EQP;
   assert((!nonlinear_mode) || (nlin_handle != NonlinearHandling::NUM_NLNHNDL));

   std::string ordering_str = config.GetOption<std::string>("model_reduction/ordering", "domain");
   if (ordering_str == "domain")    ordering = ROMOrderBy::DOMAIN;
   else if (ordering_str == "variable")  ordering = ROMOrderBy::VARIABLE;
   else
      mfem_error("ROMHandlerBase: unknown ordering!\n");
}

void ROMHandlerBase::ParseSupremizerInput(Array<int> &num_ref_supreme, Array<int> &num_supreme)
{
   assert(separate_variable);
   assert(num_var == 2); // assume vel-pres system only for now.

   num_ref_supreme.SetSize(num_rom_comp);
   num_supreme.SetSize(numSub);

   YAML::Node basis_list = config.FindNode("basis/tags");
   const int nsup_default = config.GetOption<int>("basis/number_of_supremizer", -1);

   BasisTag basis_tag;
   for (int b = 0; b < num_rom_comp; b++)
   {
      basis_tag = GetBasisTagForComponent(b, topol_handler, fom_var_names[1]);

      // By default the same number of pressure basis.
      const int nsup0 = (nsup_default > 0) ? nsup_default : num_ref_basis[b * num_var + 1];

      num_ref_supreme[b] = config.LookUpFromDict("name", basis_tag.print(), "number_of_supremizer", nsup0, basis_list);
   }

   for (int m = 0; m < numSub; m++)
   {
      int c = topol_handler->GetMeshType(m);
      num_supreme[m] = num_ref_supreme[c];
   }
}

const int ROMHandlerBase::GetBlockIndex(const int m, const int v)
{
   if (v < 0)
   {
      assert(separate_variable);
      return m;
   }

   if (!separate_variable)
      return m;
   
   if (ordering == ROMOrderBy::DOMAIN)
      return v + m * num_var;
   else if (ordering == ROMOrderBy::VARIABLE)
      return m + v * numSub;

   return -1;
}

void ROMHandlerBase::GetDomainAndVariableIndex(const int &rom_block_index, int &m, int &v)
{
   if (!separate_variable)
   {
      m = rom_block_index;
      v = -1;
      return;
   }

   if (ordering == ROMOrderBy::DOMAIN)
   {
      m = rom_block_index / num_var;
      v = rom_block_index % num_var;
   }
   else if (ordering == ROMOrderBy::VARIABLE)
   {
      m = rom_block_index % numSub;
      v = rom_block_index / numSub;
   }
}

void ROMHandlerBase::LoadReducedBasis()
{
   if (basis_loaded) return;

   std::string basis_name = basis_prefix + "_";
   int numRowRB, numColumnRB;

   carom_ref_basis.SetSize(num_rom_ref_blocks);
   basis_tags.resize(num_rom_ref_blocks);
   for (int k = 0; k < num_rom_ref_blocks; k++)
   {
      if (separate_variable)
         basis_tags[k] = GetBasisTagForComponent(k / num_var, topol_handler, fom_var_names[k % num_var]);
      else
         basis_tags[k] = GetBasisTagForComponent(k, topol_handler);

      /*
         TODO(kevin): this is a boilerplate for parallel POD/EQP training.
         We load the basis in parallel, and gather at all processes again.
         Once fully parallelized, each process will load only local part of the basis.
      */
      {
         int local_dim = dim_ref_basis[k];
         basis_reader = new CAROM::BasisReader(basis_name + basis_tags[k].print() + ".000000", CAROM::Database::formats::HDF5, local_dim, MPI_COMM_NULL);

         carom_ref_basis[k] = new CAROM::Matrix(*basis_reader->getSpatialBasis(num_ref_basis[k]));
      }
      numRowRB = carom_ref_basis[k]->numRows();
      numColumnRB = carom_ref_basis[k]->numColumns();
      printf("spatial basis-%d dimension is %d x %d\n", k, numRowRB, numColumnRB);

      delete basis_reader;
   }

   basis_loaded = true;
}

int ROMHandlerBase::GetRefIndexForSubdomain(const int &subdomain_index)
{
   int idx = topol_handler->GetMeshType(subdomain_index);
   assert(idx >= 0);
   return idx;
}

void ROMHandlerBase::SetBlockSizes()
{
   rom_block_offsets.SetSize(num_rom_blocks+1);
   rom_comp_block_offsets.SetSize(num_rom_ref_blocks+1);
   
   rom_block_offsets = 0;
   rom_comp_block_offsets = 0;

   if (separate_variable)
   {
      for (int c = 0, idx = 0; c < num_rom_comp; c++)
         for (int v = 0; v < num_var; v++, idx++)
            rom_comp_block_offsets[idx+1] = num_ref_basis[idx];

      int idx;
      for (int k = 0; k < numSub; k++)
      {
         int c = topol_handler->GetMeshType(k);
         for (int v = 0; v < num_var; v++)
         {
            idx = GetBlockIndex(k, v);

            rom_block_offsets[idx+1] = num_ref_basis[c * num_var + v];
         }  // for (int v = 0; v < num_var; v++)
      }  // for (int k = 0; k < numSub; k++)
   }
   else
   {
      for (int c = 0; c < num_rom_comp; c++)
         rom_comp_block_offsets[c+1] = num_ref_basis[c];

      for (int k = 0; k < numSub; k++)
      {
         int c = topol_handler->GetMeshType(k);
         rom_block_offsets[k+1] = num_ref_basis[c];
      }
   }

   rom_block_offsets.GetSubArray(1, num_rom_blocks, num_basis);

   rom_block_offsets.PartialSum();
   rom_comp_block_offsets.PartialSum();
}

/*
   MFEMROMHandler
*/

MFEMROMHandler::MFEMROMHandler(
   TopologyHandler *input_topol, const Array<int> &input_var_offsets,
   const std::vector<std::string> &var_names, const bool separate_variable_basis)
   : ROMHandlerBase(input_topol, input_var_offsets, var_names, separate_variable_basis)
{
   romMat = new BlockMatrix(rom_block_offsets);
   romMat->owns_blocks = true;

   std::string solver_type_str = config.GetOption<std::string>("model_reduction/linear_solver_type", "cg");
   if (solver_type_str == "direct")       linsol_type = MFEMROMHandler::SolverType::DIRECT;
   else if (solver_type_str == "cg")      linsol_type = MFEMROMHandler::SolverType::CG;
   else if (solver_type_str == "minres")      linsol_type = MFEMROMHandler::SolverType::MINRES;
   else
   {
      mfem_error("Unknown ROM linear solver type!\n");
   }

   std::string hypre_matrix_str = config.GetOption<std::string>("model_reduction/hypre_matrix", "no");
   if (hypre_matrix_str == "yes") hypre_assemble = true;
   
   if (linsol_type == MFEMROMHandler::SolverType::DIRECT)
   {
      std::string mat_type_str = config.GetOption<std::string>("model_reduction/linear_system_type", "spd");
      if (mat_type_str == "spd")          mat_type = MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE;
      else if (mat_type_str == "sid")     mat_type = MUMPSSolver::MatType::SYMMETRIC_INDEFINITE;
      else if (mat_type_str == "us")      mat_type = MUMPSSolver::MatType::UNSYMMETRIC;
      else
      {
         mfem_error("Unknown ROM linear system type!\n");
      }

      if (mat_type == MUMPSSolver::MatType::SYMMETRIC_INDEFINITE)
         mfem_warning("MUMPS matrix type SYMMETRIC_INDEFINITE can be unstable, returning inaccurate answer.\n");
   }
}

MFEMROMHandler::~MFEMROMHandler()
{
   DeletePointers(ref_basis);
   delete romMat;
   delete romMat_mono;
   delete romMat_hypre;
   delete mumps;
}

void MFEMROMHandler::LoadReducedBasis()
{
   ROMHandlerBase::LoadReducedBasis();

   ref_basis.SetSize(carom_ref_basis.Size());
   ref_basis = NULL;
   for (int k = 0; k < ref_basis.Size(); k++)
   {
      assert(carom_ref_basis[k] != NULL);
      ref_basis[k] = new DenseMatrix(carom_ref_basis[k]->numRows(), carom_ref_basis[k]->numColumns());
      CAROM::CopyMatrix(*carom_ref_basis[k], *ref_basis[k]);
   }

   dom_basis.SetSize(num_rom_blocks);
   int m, c, v, idx;
   for (int k = 0; k < num_rom_blocks; k++)
   {
      GetDomainAndVariableIndex(k, m, v);
      
      c = topol_handler->GetMeshType(m);
      idx = (separate_variable) ? c * num_var + v : c;
      dom_basis[k] = ref_basis[idx];
   }

   basis_loaded = true;
}

void MFEMROMHandler::GetReferenceBasis(const int &basis_index, DenseMatrix* &basis)
{
   assert(num_rom_ref_blocks > 0);
   assert((basis_index >= 0) && (basis_index < num_rom_ref_blocks));

   basis = ref_basis[basis_index];
}

void MFEMROMHandler::GetDomainBasis(const int &basis_index, DenseMatrix* &basis)
{
   assert(basis_loaded);
   assert((basis_index >= 0) && (basis_index < num_rom_blocks));

   basis = dom_basis[basis_index];
}

void MFEMROMHandler::ProjectOperatorOnReducedBasis(const Array2D<Operator*> &mats)
{
   assert(mats.NumRows() == num_rom_blocks);
   assert(mats.NumCols() == num_rom_blocks);
   assert(romMat->NumRowBlocks() == num_rom_blocks);
   assert(romMat->NumColBlocks() == num_rom_blocks);
   assert(basis_loaded);

   Array2D<SparseMatrix *> rom_mats;
   ProjectGlobalToDomainBasis(mats, rom_mats);
   for (int i = 0; i < num_rom_blocks; i++)
      for (int j = 0; j < num_rom_blocks; j++)
         romMat->SetBlock(i, j, rom_mats(i, j));

   romMat->Finalize();
   SetRomMat(romMat);
}

void MFEMROMHandler::SaveOperator(const std::string filename)
{
   assert(save_operator == ROMBuildingLevel::GLOBAL);
   assert(operator_loaded && romMat);

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   assert(file_id >= 0);

   hdf5_utils::WriteBlockMatrix(file_id, "ROM_matrix", romMat);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
}

void MFEMROMHandler::ProjectToRefBasis(const int &i, const Vector &vec, Vector &rom_vec)
{
   assert(basis_loaded);
   assert((i >= 0) && (i < num_rom_ref_blocks));
   DenseMatrix* basis_i;

   GetReferenceBasis(i, basis_i);
   assert(vec.Size() == basis_i->NumRows());
   assert(rom_vec.Size() == basis_i->NumCols());

   basis_i->MultTranspose(vec, rom_vec);
}

void MFEMROMHandler::ProjectToDomainBasis(const int &i, const Vector &vec, Vector &rom_vec)
{
   assert(basis_loaded);
   assert((i >= 0) && (i < num_rom_blocks));
   DenseMatrix* basis_i;

   GetDomainBasis(i, basis_i);
   assert(vec.Size() == basis_i->NumRows());
   assert(rom_vec.Size() == basis_i->NumCols());

   basis_i->MultTranspose(vec, rom_vec);
}

void MFEMROMHandler::ProjectGlobalToDomainBasis(const BlockVector* vec, BlockVector*& rom_vec)
{
   if (separate_variable)
     assert(vec->NumBlocks() == num_rom_blocks);
   else
     assert(vec->NumBlocks() == num_rom_blocks_local);

   // reset rom_vec if initiated a priori.
   if (rom_vec) delete rom_vec;

   // TODO: distribute this!
   rom_vec = new BlockVector(rom_block_offsets);

   int m, v, fom_idx;
   for (int i = localBlocks[0]; i < localBlocks[1]; ++i)
   {
      GetDomainAndVariableIndex(i, m, v);
      const int mloc = topol_handler->LocalSubdomainIndex(m);
      fom_idx = (separate_variable)? v + mloc * num_var : mloc;
      
      ProjectToDomainBasis(i, vec->GetBlock(fom_idx), rom_vec->GetBlock(i));
   }
}

void MFEMROMHandler::LiftUpFromRefBasis(const int &i, const Vector &rom_vec, Vector &vec)
{
   assert(basis_loaded);
   assert((i >= 0) && (i < num_rom_ref_blocks));
   DenseMatrix* basis_i;

   GetReferenceBasis(i, basis_i);
   assert(vec.Size() == basis_i->NumRows());
   assert(rom_vec.Size() == basis_i->NumCols());

   basis_i->Mult(rom_vec, vec);
}

void MFEMROMHandler::LiftUpFromDomainBasis(const int &i, const Vector &rom_vec, Vector &vec)
{
   assert(basis_loaded);
   assert((i >= 0) && (i < num_rom_blocks));
   DenseMatrix* basis_i;

   GetDomainBasis(i, basis_i);
   assert(vec.Size() == basis_i->NumRows());
   assert(rom_vec.Size() == basis_i->NumCols());

   basis_i->Mult(rom_vec, vec);
}

void MFEMROMHandler::LiftUpGlobal(const BlockVector &rom_vec, BlockVector &vec)
{
   assert(rom_vec.NumBlocks() == num_rom_blocks);
   if (separate_variable)
     assert(vec.NumBlocks() == num_rom_blocks);
   else
     assert(vec.NumBlocks() == num_rom_blocks_local);

   int m, v, fom_idx;
   for (int i = localBlocks[0]; i < localBlocks[1]; ++i)
   {
      GetDomainAndVariableIndex(i, m, v);
      const int mloc = topol_handler->LocalSubdomainIndex(m);
      fom_idx = (separate_variable)? v + mloc * num_var : mloc;

      LiftUpFromDomainBasis(i, rom_vec.GetBlock(i), vec.GetBlock(fom_idx));
   }
}

void MFEMROMHandler::Solve(Vector &rhs, Vector &sol)
{
   assert(operator_loaded);

   int maxIter = config.GetOption<int>("solver/max_iter", 10000);
   double rtol = config.GetOption<double>("solver/relative_tolerance", 1.e-15);
   double atol = config.GetOption<double>("solver/absolute_tolerance", 1.e-15);
   int print_level = config.GetOption<int>("solver/print_level", 0);   
   std::string prec_str = config.GetOption<std::string>("model_reduction/preconditioner", "none");

   if (linsol_type == SolverType::DIRECT)
   {
      if (topol_handler->GetRank() == 0)
	cout << "MFEMROMHandler::Solve direct, size " << rhs.Size() << endl;
      assert(mumps);
      mumps->SetPrintLevel(print_level);
      mumps->Mult(rhs, sol);
   }
   else
   {
      IterativeSolver *solver = SetIterativeSolver(linsol_type, prec_str);
      HypreParMatrix *parRomMat = NULL;
      Solver *M = NULL;    // preconditioner.
      Operator *K = NULL;  // operator.
      HypreBoomerAMG *amgM = NULL;
      // GSSmoother *gsM = NULL;

      if (prec_str == "amg")
      {
         // TODO: need to change when the actual parallelization is implemented.
         HYPRE_BigInt glob_size = rom_block_offsets.Last();
         HYPRE_BigInt row_starts[2] = {0, rom_block_offsets.Last()};
         parRomMat = new HypreParMatrix(MPI_COMM_WORLD, glob_size, row_starts, romMat_mono);
         K = parRomMat;
      }
      else if ((prec_str == "gs") || (prec_str == "none"))
         K = romMat_mono;
      else
         K = romMat;

      if (prec_str == "amg")
      {
         amgM = new HypreBoomerAMG(*parRomMat);
         amgM->SetPrintLevel(print_level);
         M = amgM;
      }
      else if (prec_str == "gs")
      {
         M = new GSSmoother(*romMat_mono);
      }
      else if (prec_str == "block_gs")
      {
         M = new BlockGSSmoother(*romMat);
      }
      else if (prec_str == "block_jacobi")
      {
         M = new BlockDSmoother(*romMat);
      }
      else if (prec_str != "none")
      {
         mfem_error("Unknown preconditioner for ROM!\n");
      }

      if (prec_str != "none")
         solver->SetPreconditioner(*M);

      if (romMat_hypre)
	solver->SetOperator(*romMat_hypre);
      else
	solver->SetOperator(*K);

      solver->SetAbsTol(atol);
      solver->SetRelTol(rtol);
      solver->SetMaxIter(maxIter);
      solver->SetPrintLevel(print_level);

      // StopWatch solveTimer;
      // solveTimer.Start();
      solver->Mult(rhs, sol);
      // solveTimer.Stop();
      // printf("ROM-solve-only time: %f seconds.\n", solveTimer.RealTime());

      // delete the created objects.
      if (prec_str == "amg")
         delete parRomMat;
      delete M;
      delete solver;
   }
}

void MFEMROMHandler::Solve(BlockVector* U)
{
   // TODO: reduced_rhs and reduced_sol are still global size!
   if (separate_variable)
     assert(U->NumBlocks() == num_rom_blocks);
   else
     assert(U->NumBlocks() == num_rom_blocks_local);
   assert(reduced_rhs);

   printf("Solve ROM.\n");
   reduced_sol = new BlockVector(rom_block_offsets);  // TODO: distribute this in parallel.
   (*reduced_sol) = 0.0;

   if (reduced_rhs_hypre.Size() > 0)
     {
       double profile_start = MPI_Wtime();
       Vector reduced_sol_hypre(reduced_rhs_hypre.Size());

       // TODO: eliminate global RHS vector in parallel.
       for (int i=0; i<reduced_rhs_hypre.Size(); ++i)
	 {
	   reduced_rhs_hypre[i] = (*reduced_rhs)[hypre_start + i];
	 }
       PrintParallelProfile("rom.solve_extract_local_rhs",
                            MPI_Wtime() - profile_start,
                            reduced_rhs_hypre.Size());

       profile_start = MPI_Wtime();
       Solve(reduced_rhs_hypre, reduced_sol_hypre);
       PrintParallelProfile("rom.solve_mumps", MPI_Wtime() - profile_start,
                            reduced_rhs_hypre.Size());

       // Gather local solutions into global solution vector.
       // TODO: keep this distributed in the parallel case.
       profile_start = MPI_Wtime();
       CAROM::Vector globalSol(reduced_sol_hypre.GetData(), reduced_sol_hypre.Size(), true);
       globalSol.gather();

       MFEM_VERIFY(globalSol.dim() == reduced_sol->Size(), "");

       for (int i=0; i<reduced_sol->Size(); ++i)
	 (*reduced_sol)[i] = globalSol(i);
       PrintParallelProfile("rom.solve_gather_global_solution",
                            MPI_Wtime() - profile_start, reduced_sol->Size());
     }
   else
     Solve(*reduced_rhs, *reduced_sol);

   // 23. reconstruct FOM state
   // TODO: distribute this in the parallel case.
   double profile_start = MPI_Wtime();
   LiftUpGlobal(*reduced_sol, *U);
   PrintParallelProfile("rom.solve_lift_solution",
                        MPI_Wtime() - profile_start, num_rom_blocks_local);
}

void MFEMROMHandler::NonlinearSolve(Operator &oper, BlockVector* U, Solver *prec)
{
   assert(U->NumBlocks() == num_rom_blocks);

   printf("Solve ROM.\n");
   reduced_sol = new BlockVector(rom_block_offsets);
   bool use_restart = config.GetOption<bool>("solver/use_restart", false);
   double amp = config.GetOption<double>("solver/initial_random_perturbation", 1.0e-5);
   if (use_restart)
      ProjectGlobalToDomainBasis(U, reduced_sol);
   else
   {
      for (int k = 0; k < reduced_sol->Size(); k++)
         (*reduced_sol)(k) = amp * UniformRandom();
   }

   int maxIter = config.GetOption<int>("solver/max_iter", 100);
   double rtol = config.GetOption<double>("solver/relative_tolerance", 1.e-10);
   double atol = config.GetOption<double>("solver/absolute_tolerance", 1.e-10);
   int print_level = config.GetOption<int>("solver/print_level", 0);   

   int jac_maxIter = config.GetOption<int>("solver/jacobian/max_iter", 10000);
   double jac_rtol = config.GetOption<double>("solver/jacobian/relative_tolerance", 1.e-10);
   double jac_atol = config.GetOption<double>("solver/jacobian/absolute_tolerance", 1.e-10);
   int jac_print_level = config.GetOption<int>("solver/jacobian/print_level", -1);
   std::string prec_str = config.GetOption<std::string>("model_reduction/preconditioner", "none");
   if (prec_str != "none") assert(prec);

   Solver *J_solver = NULL;
   if (linsol_type == SolverType::DIRECT)
   {
      mumps = new MUMPSSolver(MPI_COMM_WORLD);
      mumps->SetMatrixSymType(mat_type);
      mumps->SetPrintLevel(jac_print_level);
      J_solver = mumps;
   }
   else
   {
      IterativeSolver *iter_solver = SetIterativeSolver(linsol_type, prec_str);
      iter_solver->SetAbsTol(jac_atol);
      iter_solver->SetRelTol(jac_rtol);
      iter_solver->SetMaxIter(jac_maxIter);
      iter_solver->SetPrintLevel(jac_print_level);
      if (prec) iter_solver->SetPreconditioner(*prec);
      J_solver = iter_solver;
   }

   std::string nlin_solver = config.GetOption<std::string>("model_reduction/nonlinear_solver_type", "newton");

   if (nlin_solver == "newton")
   {
      NewtonSolver newton_solver;
      newton_solver.SetSolver(*J_solver);
      newton_solver.SetOperator(oper);
      newton_solver.SetPrintLevel(print_level); // print Newton iterations
      newton_solver.SetRelTol(rtol);
      newton_solver.SetAbsTol(atol);
      newton_solver.SetMaxIter(maxIter);

      newton_solver.Mult(*reduced_rhs, *reduced_sol);
   }
   else if (nlin_solver == "cg")
   {
      CGOptimizer optim;
      optim.SetOperator(oper);
      optim.SetPrintLevel(print_level); // print Newton iterations
      optim.SetRelTol(rtol);
      optim.SetAbsTol(atol);
      optim.SetMaxIter(maxIter);

      optim.Mult(*reduced_rhs, *reduced_sol);
   }
   else
      mfem_error("MFEMROMHandler::NonlinearSolve- Unknown ROM nonlinear solver type!\n");

   LiftUpGlobal(*reduced_sol, *U);
}

SparseMatrix* MFEMROMHandler::ProjectToRefBasis(const int &i, const int &j, const Operator *mat)
{
   assert((i >= 0) && (i < num_rom_ref_blocks));
   assert((j >= 0) && (j < num_rom_ref_blocks));
   // assert(mat->Finalized());
   assert(basis_loaded);
   
   DenseMatrix *basis_i, *basis_j;
   GetReferenceBasis(i, basis_i);
   GetReferenceBasis(j, basis_j);

   assert(basis_i->NumRows() == mat->Height());
   assert(basis_j->NumRows() == mat->Width());

   if (mat)
      return mfem::SparseRtAP(*basis_i, *mat, *basis_j);
   else
      return new SparseMatrix(basis_i->NumCols(), basis_j->NumCols());
}

SparseMatrix* MFEMROMHandler::ProjectToDomainBasis(const int &i, const int &j, const Operator *mat)
{
   assert((i >= 0) && (i < num_rom_blocks));
   assert((j >= 0) && (j < num_rom_blocks));
   // assert(mat->Finalized());
   assert(basis_loaded);
   
   DenseMatrix *basis_i, *basis_j;
   GetDomainBasis(i, basis_i);
   GetDomainBasis(j, basis_j);

   if (mat)
   {
      assert(basis_i->NumRows() == mat->Height());
      assert(basis_j->NumRows() == mat->Width());
      return mfem::SparseRtAP(*basis_i, *mat, *basis_j);
   }
   else 
      return new SparseMatrix(basis_i->NumCols(), basis_j->NumCols());
}

void MFEMROMHandler::ProjectToRefBasis(
   const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   assert(basis_loaded);
   assert((idx_i.Min() >= 0) && (idx_i.Max() < num_rom_ref_blocks));
   assert((idx_j.Min() >= 0) && (idx_j.Max() < num_rom_ref_blocks));
   assert(idx_i.Size() == mats.NumRows());
   assert(idx_j.Size() == mats.NumCols());

   /* Clean up existing rom_mats */
   for (int i = 0; i < rom_mats.NumRows(); i++)
      for (int j = 0; j < rom_mats.NumCols(); j++)
         if (rom_mats(i, j)) delete rom_mats(i, j);
   
   rom_mats.SetSize(mats.NumRows(), mats.NumCols());
   for (int i = 0; i < rom_mats.NumRows(); i++)
      for (int j = 0; j < rom_mats.NumCols(); j++)
         rom_mats(i, j) = ProjectToRefBasis(idx_i[i], idx_j[j], mats(i, j));
}

void MFEMROMHandler::ProjectToDomainBasis(
   const Array<int> &idx_i, const Array<int> &idx_j, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   assert(basis_loaded);
   assert((idx_i.Min() >= 0) && (idx_i.Max() < num_rom_blocks));
   assert((idx_j.Min() >= 0) && (idx_j.Max() < num_rom_blocks));
   assert(idx_i.Size() == mats.NumRows());
   assert(idx_j.Size() == mats.NumCols());

   /* Clean up existing rom_mats */
   for (int i = 0; i < rom_mats.NumRows(); i++)
      for (int j = 0; j < rom_mats.NumCols(); j++)
         if (rom_mats(i, j)) delete rom_mats(i, j);
   
   rom_mats.SetSize(mats.NumRows(), mats.NumCols());
   for (int i = 0; i < rom_mats.NumRows(); i++)
      for (int j = 0; j < rom_mats.NumCols(); j++)
         rom_mats(i, j) = ProjectToDomainBasis(idx_i[i], idx_j[j], mats(i, j));
}

void MFEMROMHandler::ProjectComponentToRefBasis(
   const int &c, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   int size = (separate_variable) ? num_var : 1;
   Array<int> idx(size);
   idx = c;
   if (separate_variable)
      for (int v = 0; v < num_var; v++)
      {
         idx[v] *= num_var;
         idx[v] += v;
      }

   assert((mats.NumRows() == size) && (mats.NumCols() == size));

   ProjectToRefBasis(idx, idx, mats, rom_mats);
}

void MFEMROMHandler::ProjectComponentToDomainBasis(
   const int &c, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   int size = (separate_variable) ? num_var : 1;
   Array<int> idx(size);
   idx = c;
   if (separate_variable)
      for (int v = 0; v < num_var; v++)
      {
         idx[v] *= num_var;
         idx[v] += v;
      }

   assert((mats.NumRows() == size) && (mats.NumCols() == size));

   ProjectToDomainBasis(idx, idx, mats, rom_mats);
}

void MFEMROMHandler::ProjectInterfaceToRefBasis(
   const int &c1, const int &c2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   int size = (separate_variable) ? 2 * num_var : 2;
   Array<int> idx(size);
   if (separate_variable)
      for (int v = 0; v < num_var; v++)
      {
         idx[v] = c1 * num_var + v;
         idx[v + num_var] = c2 * num_var + v;
      }
   else
   {
      idx[0] = c1;
      idx[1] = c2;
   }

   assert((mats.NumRows() == size) && (mats.NumCols() == size));

   ProjectToRefBasis(idx, idx, mats, rom_mats);
}

void MFEMROMHandler::ProjectInterfaceToDomainBasis(
   const int &c1, const int &c2, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   int size = (separate_variable) ? 2 * num_var : 2;
   Array<int> idx(size);
   if (separate_variable)
      for (int v = 0; v < num_var; v++)
      {
         idx[v] = c1 * num_var + v;
         idx[v + num_var] = c2 * num_var + v;
      }
   else
   {
      idx[0] = c1;
      idx[1] = c2;
   }

   assert((mats.NumRows() == size) && (mats.NumCols() == size));

   ProjectToDomainBasis(idx, idx, mats, rom_mats);
}

void MFEMROMHandler::ProjectVariableToDomainBasis(
   const int &vi, const int &vj, const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   assert(separate_variable);

   int size = numSub;
   Array<int> idx_i(size), idx_j(size);
   for (int m = 0; m < numSub; m++)
   {
      idx_i[m] = m * num_var + vi;
      idx_j[m] = m * num_var + vj;
   }

   assert((mats.NumRows() == size) && (mats.NumCols() == size));

   ProjectToDomainBasis(idx_i, idx_j, mats, rom_mats);
}

void MFEMROMHandler::ProjectGlobalToDomainBasis(const Array2D<Operator*> &mats, Array2D<SparseMatrix *> &rom_mats)
{
   assert((mats.NumRows() == num_rom_blocks) && (mats.NumCols() == num_rom_blocks));

   Array<int> idx(num_rom_blocks);
   for (int k = 0; k < num_rom_blocks; k++) idx[k] = k;

   ProjectToDomainBasis(idx, idx, mats, rom_mats);
}

void MFEMROMHandler::LoadOperatorFromFile(const std::string filename)
{
   assert(save_operator == ROMBuildingLevel::GLOBAL);
   delete romMat;
   delete romMat_mono;

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);   

   romMat = hdf5_utils::ReadBlockMatrix(file_id, "ROM_matrix", rom_block_offsets);
   SetRomMat(romMat);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
}

// TODO: should this be in the base class ROMHandlerBase?
void MFEMROMHandler::LoadBalanceROMBlocks(int rank, int nproc)
{
  MFEM_VERIFY(num_rom_blocks + 1 == rom_block_offsets.Size(), "");
  MFEM_VERIFY(num_rom_blocks >= nproc, "");

  const int gsize = rom_block_offsets[num_rom_blocks];

  const int bp = num_rom_blocks / nproc;  // Number of blocks per rank
  const int ne = num_rom_blocks - (bp * nproc);  // Number of ranks needing an extra block

  localSizes.SetSize(nproc);
  localSizes = 0;

  localNumBlocks.SetSize(nproc);
  localNumBlocks = 0;

  localBlocks[0] = 0;

  int sum = 0;
  int bsum = 0;
  for (int j=0, i=0; j<nproc; ++j)
    {
      const int nb = j < ne ? bp + 1 : bp;  // Number of blocks for rank j
      localNumBlocks[j] = nb;
      for (int b=0; b<nb; ++b)
	{
	  const int bsize = rom_block_offsets[i + 1] - rom_block_offsets[i];
	  localSizes[j] += bsize;
	  sum += bsize;
	  i++;
	}

      bsum += nb;

      if (j == rank)
	localBlocks[1] = bsum;
      else if (j == rank - 1)
	localBlocks[0] = bsum;
    }

  MFEM_VERIFY(sum == gsize && bsum == num_rom_blocks, "");

}

void MFEMROMHandler::CreateHypreParMatrix(BlockMatrix *input_mat, int rank, int nproc)
{
  double profile_start = MPI_Wtime();
  const HYPRE_BigInt gsize = rom_block_offsets[num_rom_blocks];

  Array<int> boffset(nproc+1);

  Array<HYPRE_BigInt> offset(nproc+1);
  boffset[0] = 0;
  offset[0] = 0;
  for (int i=1; i<=nproc; ++i)
    {
      offset[i] = offset[i - 1] + localSizes[i - 1];
      boffset[i] = boffset[i - 1] + localNumBlocks[i - 1];
    }

  const int loc0 = offset[rank];
  const int loc1 = offset[rank] + localSizes[rank];

  std::set<int> offd_blocks;
  std::set<int> offd_cols;

  const int nblocks = localNumBlocks[rank];
  num_rom_blocks_local = nblocks;

  for (int b=0; b<nblocks; ++b)
    {
      // TODO: eliminate global blocks (just assemble local blocks).
      const int gb = boffset[rank] + b; // Global block index
      for (int j=0; j<num_rom_blocks; ++j)
	{
	  if (boffset[rank] <= j && j < boffset[rank] + nblocks)
	    continue;

	  // TODO: store sparse blocks, instead of looping over all blocks to check for nonzeros
	  if (!input_mat->IsZeroBlock(gb, j))
	    {
	      offd_blocks.insert(j);
	    }
	}
    }
  PrintParallelProfile("rom.hypre_discover_offd_blocks",
                       MPI_Wtime() - profile_start,
                       static_cast<long>(nblocks) * num_rom_blocks);

  int offd_size = 0;
  for (auto b : offd_blocks)
    {
      const int bsize_j = rom_block_offsets[b + 1] - rom_block_offsets[b];
      offd_size += bsize_j;
    }

  // The following construction of a HypreParMatrix, by using a SparseMatrix
  // for diag and offd, follows the example of
  // ParFiniteElementSpace::ParallelDerefinementMatrix
  hdiag = new SparseMatrix(localSizes[rank], localSizes[rank]);
  SparseMatrix &diag = *hdiag;
  hoffd = new SparseMatrix(localSizes[rank], gsize);
  SparseMatrix &offd = *hoffd;

  // Set diag and offd

  profile_start = MPI_Wtime();
  int localOffset = 0;
  for (int b=0; b<nblocks; ++b)
    {
      // TODO: eliminate global blocks (just assemble local blocks).
      const int gb = boffset[rank] + b; // Global block index

      const int bsize = rom_block_offsets[gb + 1] - rom_block_offsets[gb];
      Array<int> rows(bsize);
      for (int i=0; i<bsize; ++i)
	{
	  rows[i] = localOffset + i;
	}

      localOffset += bsize;

      for (int j=0; j<num_rom_blocks; ++j)
	{
	  // TODO: store sparse blocks, instead of looping over all blocks to check for nonzeros
	  if (!input_mat->IsZeroBlock(gb, j))
	    {
	      const SparseMatrix &block = input_mat->GetBlock(gb, j);
	      // TODO: this conversion to DenseMatrix is inefficient. If ROM blocks are always
	      // dense, can we just store them as DenseMatrix instances in the first place?

	      DenseMatrix *db = block.ToDenseMatrix();

	      const bool diagBlock = boffset[rank] <= j && j < boffset[rank] + nblocks;

	      const int bsize_j = rom_block_offsets[j + 1] - rom_block_offsets[j];
	      Array<int> cols(bsize_j);
	      for (int i=0; i<bsize_j; ++i)
		{
		  cols[i] = rom_block_offsets[j] + i;

		  // TODO: simplify
		  MFEM_VERIFY((cols[i] < loc0 || cols[i] >= loc1) == !diagBlock, "");

		  if (cols[i] < loc0 || cols[i] >= loc1)  // Off-diagonal
		    {
		      offd_cols.insert(cols[i]);
		    }

		  if (diagBlock)
		    cols[i] -= loc0;
		}

	      if (diagBlock)
		{
		  // Diagonal block
		  MFEM_VERIFY(bsize_j == bsize, "");
		  diag.AddSubMatrix(rows, cols, *db);
		}
	      else
		{
		  // Off-diagonal block
		  offd.AddSubMatrix(rows, cols, *db);
		}

	      delete db;
	    }
	}
    }
  PrintParallelProfile("rom.hypre_fill_diag_offd",
                       MPI_Wtime() - profile_start,
                       static_cast<long>(nblocks) * num_rom_blocks);

  const int num_offd_cols = offd_cols.size();
  MFEM_VERIFY(num_offd_cols == offd_size, "");
  cmap = new HYPRE_BigInt[num_offd_cols];
  std::map<int, int> cmap_inv;

  int cnt = 0;
  for (auto col : offd_cols)
    {
      cmap[cnt] = col;
      cmap_inv[col] = cnt;

      cnt++;
    }

  MFEM_VERIFY(cnt == num_offd_cols, "");

  profile_start = MPI_Wtime();
  diag.Finalize();
  offd.Finalize();

  if (num_offd_cols > 0)
    {
      // Map column indices in offd
      int *offI = offd.GetI();
      int *offJ = offd.GetJ();

      const int ne = offI[localSizes[rank]];  // Total number of entries in offd

      for (int i=0; i<ne; ++i)
	{
	  const int c = cmap_inv[offJ[i]];
	  offJ[i] = c;
	}

      offd.SortColumnIndices();
    }

  offd.SetWidth(offd_size);
  PrintParallelProfile("rom.hypre_finalize_and_remap",
                       MPI_Wtime() - profile_start, num_offd_cols);

  Array<HYPRE_BigInt> starts(2);
  starts[0] = offset[rank];
  starts[1] = offset[rank + 1];

  MFEM_VERIFY(HYPRE_AssumedPartitionCheck(), "");

  profile_start = MPI_Wtime();
  if (nproc == 1) // Serial case
    {
      // constructor with 4 arguments, v1
      romMat_hypre = new HypreParMatrix(MPI_COMM_WORLD, gsize,
					starts.GetData(), &diag);
    }
  else
    {
      // constructor with 8+1 arguments
      romMat_hypre = new HypreParMatrix(MPI_COMM_WORLD, gsize, gsize,
					starts.GetData(), starts.GetData(),
					hdiag, hoffd, cmap, true);
      romMat_hypre->SetOwnerFlags(romMat_hypre->OwnsDiag(), romMat_hypre->OwnsOffd(), 1);
    }
  PrintParallelProfile("rom.hypre_matrix_constructor",
                       MPI_Wtime() - profile_start, localSizes[rank]);

  hypre_start = starts[0];

  reduced_rhs_hypre.SetSize(starts[1] - starts[0]);

  profile_start = MPI_Wtime();
  delete mumps;
  mumps = new MUMPSSolver(MPI_COMM_WORLD);
  mumps->SetMatrixSymType(mat_type);
  mumps->SetOperator(*romMat_hypre);
  PrintParallelProfile("rom.mumps_setup", MPI_Wtime() - profile_start,
                       localSizes[rank]);
}

void MFEMROMHandler::SetRomMat(BlockMatrix *input_mat, const bool init_direct_solver)
{
   if (romMat != input_mat)
   {
      delete romMat;
      romMat = input_mat;
   }

   delete romMat_mono;
   romMat_mono = NULL;
   if (init_direct_solver)
      romMat_mono = romMat->CreateMonolithic();

   if ((linsol_type == SolverType::DIRECT) && (init_direct_solver))
      SetupDirectSolver();
   operator_loaded = true;
}

void MFEMROMHandler::SaveRomSystem(const std::string &input_prefix, const std::string type)
{
   if (topol_handler->GetRank() != 0) return; // Only the root process writes files.

   if (!romMat_mono)
   {
      assert(romMat);
      romMat_mono = romMat->CreateMonolithic();
   }
   assert(reduced_rhs);
   assert(reduced_sol);

   PrintVector(*reduced_rhs, input_prefix + "_rhs.txt");
   PrintVector(*reduced_sol, input_prefix + "_sol.txt");

   std::string matfile = input_prefix + "_mat." + type;
   std::ofstream file(matfile.c_str());
   if (type == "mm")
      romMat_mono->PrintMM(file);
   else if (type == "matlab")
      romMat_mono->PrintMatlab(file);
   else
      mfem_error("MFEMROMHandler::SaveRomSystem - unknown matrix format type!\n");
}

void MFEMROMHandler::SaveBasisVisualization(
   const Array<FiniteElementSpace *> &fes, const std::vector<std::string> &var_names)
{
   if (!save_basis_visual) return;
   assert(basis_loaded);

   const int num_var = var_names.size();
   assert(fes.Size() == num_var * numSub);

   std::string visual_prefix = config.GetRequiredOption<std::string>("basis/visualization/prefix");

   for (int c = 0; c < num_rom_ref_blocks; c++)
   {
      std::string file_prefix = visual_prefix;
      file_prefix += "_" + std::to_string(c);

      int midx = -1;
      for (int m = 0; m < numSub; m++)
         if (topol_handler->GetMeshType(m) == c) { midx = m; break; }
      assert(midx >= 0);

      Mesh *mesh = fes[midx * num_var]->GetMesh();
      const int order = fes[midx * num_var]->FEColl()->GetOrder();
      ParaViewDataCollection coll = ParaViewDataCollection(file_prefix.c_str(), mesh);

      Array<int> var_offsets(num_ref_basis[c] * num_var + 1);
      var_offsets[0] = 0;
      for (int k = 0, vidx = 1; k < num_ref_basis[c]; k++)
         for (int v = 0, idx = midx * num_var; v < num_var; v++, idx++, vidx++)
            var_offsets[vidx] = fes[idx]->GetVSize();
      var_offsets.PartialSum();
      BlockVector basis_view(ref_basis[c]->GetData(), var_offsets);

      Array<GridFunction*> basis_gf(num_ref_basis[c] * num_var);
      basis_gf = NULL;
      for (int k = 0, idx = 0; k < num_ref_basis[c]; k++)
      {
         for (int v = 0, fidx = midx * num_var; v < num_var; v++, idx++, fidx++)
         {
            std::string field_name = var_names[v] + "_basis_" + std::to_string(k);
            basis_gf[idx] = new GridFunction(fes[fidx], basis_view.GetBlock(idx), 0);
            coll.RegisterField(field_name.c_str(), basis_gf[idx]);   
         }
      }

      coll.SetLevelsOfDetail(order);
      coll.SetHighOrderOutput(true);
      coll.SetPrecision(8);
      coll.SetOwnData(false);
      coll.Save();

      for (int k = 0; k < basis_gf.Size(); k++) delete basis_gf[k];
   }
}

IterativeSolver* MFEMROMHandler::SetIterativeSolver(const MFEMROMHandler::SolverType &linsol_type_, const std::string &prec_type)
{
   IterativeSolver *solver;
   switch (linsol_type_)
   {
      case (SolverType::CG):
      {
         if (prec_type == "amg") solver = new CGSolver(MPI_COMM_WORLD);
         else                    solver = new CGSolver(MPI_COMM_WORLD);
         break;
      }
      case (SolverType::MINRES):
      {
         if (prec_type == "amg") solver = new MINRESSolver(MPI_COMM_WORLD);
         else                    solver = new MINRESSolver();
         break;
      }
      case (SolverType::GMRES):
      {
         if (prec_type == "amg") solver = new GMRESSolver(MPI_COMM_WORLD);
         else                    solver = new GMRESSolver();
         break;
      }
      default:
      {
         mfem_error("Unknown ROM iterative linear solver type!\n");
         break;
      }
   }

   return solver;
}

void MFEMROMHandler::SetupDirectSolver()
{
   // If nonlinear mode, Jacobian will keep changing within Solve, thus no need of initial LU factorization.
   if ((linsol_type != MFEMROMHandler::SolverType::DIRECT))
      return;

   // TODO: need to change when the actual parallelization is implemented.
   sys_glob_size = romMat_mono->NumRows();
   sys_row_starts[0] = 0;
   sys_row_starts[1] = romMat_mono->NumRows();
   romMat_hypre = new HypreParMatrix(MPI_COMM_SELF, sys_glob_size, sys_row_starts, romMat_mono);

   mumps = new MUMPSSolver(MPI_COMM_SELF);
   mumps->SetMatrixSymType(mat_type);
   mumps->SetOperator(*romMat_hypre);
}

void MFEMROMHandler::AppendReferenceBasis(const int &idx, const DenseMatrix &mat)
{
   assert(basis_loaded);
   assert((idx >= 0) && (idx < num_rom_ref_blocks));
   assert(ref_basis[idx]->NumRows() == mat.NumRows());
   assert(ref_basis[idx]->NumCols() == num_ref_basis[idx]);

   DenseMatrix tmp(*ref_basis[idx]);

   int add_col = mat.NumCols();
   // NOTE(kevin): expanding a DenseMatrix does not preserve its data.
   ref_basis[idx]->SetSize(ref_basis[idx]->NumRows(), num_ref_basis[idx] + add_col);
   ref_basis[idx]->SetSubMatrix(0, 0, tmp);
   ref_basis[idx]->SetSubMatrix(0, num_ref_basis[idx], mat);

   num_ref_basis[idx] += add_col;

   // Reset the block offsets and domain num_basis.
   SetBlockSizes();

   printf("spatial basis-%d dimension updated: %d x %d\n",
          idx, ref_basis[idx]->NumRows(), ref_basis[idx]->NumCols());
}

// void MFEMROMHandler::GetBlockSparsity(
//    const SparseMatrix *mat, const Array<int> &block_offsets, Array2D<bool> &mat_zero_blocks)
// {
//    assert(mat);
//    assert(block_offsets.Size() > 0);
//    assert(mat->NumRows() == block_offsets.Last());
//    assert(mat->NumCols() == block_offsets.Last());

//    const int num_blocks = block_offsets.Size() - 1;
//    mat_zero_blocks.SetSize(num_blocks, num_blocks);
//    mat_zero_blocks = true;

//    DenseMatrix tmp;
//    Array<int> rows, cols;
//    for (int i = 0; i < num_blocks; i++)
//    {
//       rows.SetSize(block_offsets[i+1] - block_offsets[i]);
//       for (int ii = block_offsets[i], idx = 0; ii < block_offsets[i+1]; ii++, idx++)
//          rows[idx] = ii;

//       for (int j = 0; j < num_blocks; j++)
//       {
//          cols.SetSize(block_offsets[j+1] - block_offsets[j]);
//          for (int jj = block_offsets[j], idx = 0; jj < block_offsets[j+1]; jj++, idx++)
//             cols[idx] = jj;

//          tmp.SetSize(rows.Size(), cols.Size());
//          tmp = 0.0;
//          mat->GetSubMatrix(rows, cols, tmp);
//          mat_zero_blocks(i, j) = CheckZeroBlock(tmp);
//       }  // for (int j = 0; j < num_blocks; j++)
//    }  // for (int i = 0; i < num_blocks; i++)
// }

// bool MFEMROMHandler::CheckZeroBlock(const DenseMatrix &mat)
// {
//    for (int i = 0; i < mat.NumRows(); i++)
//       for (int j = 0; j < mat.NumCols(); j++)
//          if (mat(i, j) != 0.0)
//             return false;

//    return true;
// }

}  // namespace mfem
