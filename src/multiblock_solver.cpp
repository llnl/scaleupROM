// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "multiblock_solver.hpp"
#include "linalg_utils.hpp"
#include "component_topology_handler.hpp"
#include "etc.hpp"

using namespace std;
using namespace mfem;

MultiBlockSolver::MultiBlockSolver()
{
   /*
      TODO(kevin): this is a boilerplate for parallel POD/EQP training.
      In full parallelization, MultiBlockSolver will take MPI_Comm as an input argument.
   */
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   ParseInputs();

   TopologyData topol_data;
   switch (topol_mode)
   {
      case TopologyHandlerMode::SUBMESH:
      {
         topol_handler = new SubMeshTopologyHandler();
         break;
      }
      case TopologyHandlerMode::COMPONENT:
      {
         topol_handler = new ComponentTopologyHandler();
         break;
      }
      default:
      {
         mfem_error("Unknown topology handler mode!\n");
         break;
      }
   }
   topol_handler->ExportInfo(meshes, topol_data);
   
   // Receive topology info
   numSub = topol_data.numSub;
   numSubLoc = topol_handler->GetNumLocalSubdomains();
   numSubStored = meshes.Size();
   dim = topol_data.dim;
   global_bdr_attributes = *(topol_data.global_bdr_attributes);
   numBdr = global_bdr_attributes.Size();
}

MultiBlockSolver::~MultiBlockSolver()
{
   delete U;
   delete RHS;

   if (visual.save)
      for (int k = 0; k < paraviewColls.Size(); k++) delete paraviewColls[k];

   for (int k = 0; k < us.Size(); k++) delete us[k];
   for (int k = 0; k < fes.Size(); k++) delete fes[k];
   for (int k = 0; k < fec.Size(); k++) delete fec[k];

   for (int k = 0; k < global_fes.Size(); k++) delete global_fes[k];
   for (int k = 0; k < global_us_visual.Size(); k++) delete global_us_visual[k];
   DeletePointers(error_visual);
   DeletePointers(global_error_visual);

   for (int k = 0; k < bdr_markers.Size(); k++)
      delete bdr_markers[k];

   delete rom_handler;
   delete topol_handler;

   DeletePointers(comp_fes);
   delete rom_elems;
}

void MultiBlockSolver::ParseInputs()
{
   topol_mode = SetTopologyHandlerMode();

   order = config.GetOption<int>("discretization/order", 1);
   full_dg = config.GetOption<bool>("discretization/full-discrete-galerkin", false);

   // solver option;
   use_amg = config.GetOption<bool>("solver/use_amg", true);
   direct_solve = config.GetOption<bool>("solver/direct_solve", false);

   visual.save = config.GetOption<bool>("visualization/enabled", false);
   if (visual.save)
   {
      // Default file path if no input file name is provided.
      visual.dir = config.GetOption<std::string>("visualization/file_path/directory", ".");
      visual.prefix = config.GetOption<std::string>("visualization/file_path/prefix", "paraview_output");

      visual.unified_view = config.GetOption<bool>("visualization/unified_paraview", false);
      if (visual.unified_view)
      {
         if (topol_mode == TopologyHandlerMode::COMPONENT)
            mfem_error("ComponentTopologyHandler does not yet support unified paraview!\n");
      }
      else
      {
         visual.domain_offset = config.GetOption<int>("visualization/domain_offset", 0);
         visual.domain_interval = config.GetOption<int>("visualization/domain_frequency", 1);
      }

      // visualization type for pointwise error.
      visual.save_error = config.GetOption<bool>("visualization/visualize_error", false);

      visual.time_interval = config.GetOption<int>("visualization/time_interval", -1);
   }

   // rom inputs.
   use_rom = config.GetOption<bool>("main/use_rom", false);
   separate_variable_basis = config.GetOption<bool>("model_reduction/separate_variable_basis", false);

   // save solution if single run.
   SetSolutionSaveMode(config.GetOption<bool>("save_solution/enabled", false));
}

void MultiBlockSolver::SetSolutionSaveMode(const bool save_sol_)
{
   // save solution if single run.
   save_sol = save_sol_;
   if (save_sol)
   {
      // Default file path if no input file name is provided.
      sol_dir = config.GetOption<std::string>("save_solution/file_path/directory", ".");
      sol_prefix = config.GetOption<std::string>("save_solution/file_path/prefix", "solution");
   }
}

void MultiBlockSolver::GetVariableVector(const int &var_idx, BlockVector &global, BlockVector &var)
{
   assert((var_idx >= 0) && (var_idx < num_var));
   assert(global.NumBlocks() == (num_var * numSub));
   assert(var.NumBlocks() == (numSub));

   for (int m = 0; m < numSub; m++)
   {
      int g_idx = num_var * m + var_idx;
      assert(var.BlockSize(m) == global.BlockSize(g_idx));

      Vector tmp;
      var.GetBlockView(m, tmp);
      tmp = global.GetBlock(g_idx);
   }
}

void MultiBlockSolver::SetVariableVector(const int &var_idx, BlockVector &var, BlockVector &global)
{
   assert((var_idx >= 0) && (var_idx < num_var));
   assert(global.NumBlocks() == (num_var * numSub));
   assert(var.NumBlocks() == (numSub));

   for (int m = 0; m < numSub; m++)
   {
      int g_idx = num_var * m + var_idx;
      assert(var.BlockSize(m) == global.BlockSize(g_idx));

      Vector tmp;
      global.GetBlockView(g_idx, tmp);
      tmp = var.GetBlock(m);
   }
}

void MultiBlockSolver::SortBySubdomains(BlockVector &by_var, BlockVector &by_sub)
{
   assert(by_var.NumBlocks() == (num_var * numSubLoc));
   assert(by_sub.NumBlocks() == (num_var * numSubLoc));

   for (int m = 0; m < numSubLoc; m++)
      for (int v = 0; v < num_var; v++)
      {
         int by_var_idx = numSubLoc * v + m;
         int by_sub_idx = num_var * m + v;
         assert(by_var.BlockSize(by_var_idx) == by_sub.BlockSize(by_sub_idx));

         Vector tmp;
         by_sub.GetBlockView(by_sub_idx, tmp);
         tmp = by_var.GetBlock(by_var_idx);
      }
}

void MultiBlockSolver::SortByVariables(BlockVector &by_sub, BlockVector &by_var)
{
   assert(by_var.NumBlocks() == (num_var * numSubLoc));
   assert(by_sub.NumBlocks() == (num_var * numSubLoc));

   for (int m = 0; m < numSubLoc; m++)
      for (int v = 0; v < num_var; v++)
      {
         int by_var_idx = numSubLoc * v + m;
         int by_sub_idx = num_var * m + v;
         assert(by_var.BlockSize(by_var_idx) == by_sub.BlockSize(by_sub_idx));

         Vector tmp;
         by_var.GetBlockView(by_var_idx, tmp);
         tmp = by_sub.GetBlock(by_sub_idx);
      }
}

void MultiBlockSolver::SetupBCVariables()
{
   // Set up boundary markers.
   int max_bdr_attr = -1;
   for (int m = 0; m < numSubLoc; m++)
   {
      max_bdr_attr = max(max_bdr_attr, meshes[m]->bdr_attributes.Max());
   }
   {
     const int local_max = max_bdr_attr;
     MPI_Allreduce(&local_max, &max_bdr_attr, 1, MPI_INT, MPI_MAX, topol_handler->GetComm());
   }

   // TODO: technically this should be Array<Array2D<int>*> for each meshes.
   // Running with MFEM debug version will lead to error when assembling boundary integrators.
   bdr_markers.SetSize(global_bdr_attributes.Size());
   for (int k = 0; k < bdr_markers.Size(); k++) {
      int bdr_attr = global_bdr_attributes[k];
      assert((bdr_attr > 0) && (bdr_attr <= max_bdr_attr));
      bdr_markers[k] = new Array<int>(max_bdr_attr);
      (*bdr_markers[k]) = 0;
      (*bdr_markers[k])[bdr_attr-1] = 1;
   }

   bdr_type.SetSize(global_bdr_attributes.Size());
   bdr_type = BoundaryType::NUM_BDR_TYPE;
}

void MultiBlockSolver::SetBdrType(const BoundaryType type, const int &global_battr_idx)
{
   assert(bdr_type.Size() == global_bdr_attributes.Size());
   if (global_battr_idx < 0)
      bdr_type = type;
   else
   {
      assert(global_battr_idx < global_bdr_attributes.Size());
      bdr_type[global_battr_idx] = type;
   }
}

void MultiBlockSolver::GetComponentFESpaces(Array<FiniteElementSpace *> &comp_fes)
{
   assert(fec.Size() == num_var);
   assert(vdim.Size() == num_var);

   // Component domain system
   const int num_comp = topol_handler->GetNumComponents();
   comp_fes.SetSize(num_var * num_comp);
   comp_fes = NULL;

   for (int c = 0, idx = 0; c < num_comp; c++)
   {
      Mesh *comp = topol_handler->GetComponentMesh(c);
      for (int v = 0; v < num_var; v++, idx++)
         comp_fes[idx] = new FiniteElementSpace(comp, fec[v], vdim[v]);
   }
}

void MultiBlockSolver::BuildROMLinElems()
{
   assert(topol_mode == TopologyHandlerMode::COMPONENT);
   assert(rom_handler->BasisLoaded());

   BuildCompROMLinElems();

   // Boundary penalty matrices
   BuildBdrROMLinElems();

   // Port penalty matrices
   BuildItfaceROMLinElems();
}

void MultiBlockSolver::AssembleROMMat()
{
   double profile_start = MPI_Wtime();
   const Array<int> *rom_block_offsets = rom_handler->GetBlockOffsets();
   BlockMatrix *romMat = new BlockMatrix(*rom_block_offsets);
   romMat->owns_blocks = true;

   AssembleROMMat(*romMat);
   PrintParallelProfile("rom.block_assembly", MPI_Wtime() - profile_start,
                        numSubLoc);

   profile_start = MPI_Wtime();
   romMat->Finalize();
   PrintParallelProfile("rom.block_finalize", MPI_Wtime() - profile_start);
   const bool hypreAssemble = rom_handler->HypreAssemble();
   profile_start = MPI_Wtime();
   rom_handler->SetRomMat(romMat, !hypreAssemble);
   PrintParallelProfile("rom.set_matrix_monolithic",
                        MPI_Wtime() - profile_start);

   if (hypreAssemble)
   {
       profile_start = MPI_Wtime();
       if (!separate_variable_basis)
       {
          Array<int> local_num_blocks;
          rom_handler->GetLocalNumBlocks(local_num_blocks);
          Array<int> block_starts(nproc + 1);
          block_starts[0] = 0;
          for (int r = 0; r < nproc; r++)
             block_starts[r + 1] = block_starts[r] + local_num_blocks[r];

          SparseBlockMatrixLayout layout(*rom_block_offsets);
          SparseBlockMatrixDistribution distribution(layout, block_starts);
          Array<int> stored_rows(distribution.NumBlocks(rank));
          for (int i = 0; i < stored_rows.Size(); i++)
             stored_rows[i] = distribution.FirstBlock(rank) + i;

          SparseBlockMatrixBuilder builder(layout, stored_rows);
          AssembleROMMat(builder);
          SparseBlockMatrix sparse_matrix = builder.Finalize();
          rom_handler->CreateHypreParMatrix(layout, distribution, sparse_matrix);
       }
       else
       {
          // Separate-variable ownership still needs an explicit logical-to-
          // distributed block numbering. Retain the established path until
          // that mapping is available.
          rom_handler->CreateHypreParMatrix(romMat, rank, nproc);
       }
       PrintParallelProfile("rom.create_hypre_total",
                            MPI_Wtime() - profile_start);

       localBlocks = rom_handler->GetLocalBlocks(); // TODO: get this for FOM without depending on ROM?
       rom_handler->GetLocalNumBlocks(localNumBlocks); // TODO: get this for FOM without depending on ROM?

       topol_handler->GetNeighbors(neighbors);
   }
}

void MultiBlockSolver::AssembleROMMat(BlockMatrix &romMat)
{
   assert(topol_mode == TopologyHandlerMode::COMPONENT);
   assert(rom_elems);

   const int ossub = topol_handler->GlobalSubdomainIndex(0); // Offset for first local subdomain

   double profile_start = MPI_Wtime();
   // component domain matrix.
   for (int mm = 0; mm < numSubLoc; mm++)
     {
      const int m = mm + ossub;
      const int m_global = topol_handler->GlobalSubdomainIndex(mm);
      const int c_type = topol_handler->GetMeshType(m_global);
      const int num_block = (separate_variable_basis) ? num_var : 1;

      Array<int> midx(num_block);
      for (int v = 0; v < num_block; v++)
         midx[v] = rom_handler->GetBlockIndex(m, v);

      MatrixBlocks *comp_mat = rom_elems->comp[c_type];
      AddToBlockMatrix(midx, midx, *comp_mat, romMat);

      // boundary matrices of each component.
      Array<int> *bdr_c2g = topol_handler->GetBdrAttrComponentToGlobalMap(mm);

      for (int b = 0; b < bdr_c2g->Size(); b++)
      {
         int global_idx = global_bdr_attributes.Find((*bdr_c2g)[b]);
         if (global_idx < 0) continue;
         if (!BCExistsOnBdr(global_idx)) continue;

         /* we assume only Neumann condition would not add an operator. */
         if (bdr_type[global_idx] == BoundaryType::NEUMANN)
            continue;

         MatrixBlocks *bdr_mat = (*(rom_elems->bdr[c_type]))[b];
         AddToBlockMatrix(midx, midx, *bdr_mat, romMat);
      }  // for (int b = 0; b < bdr_c2g->Size(); b++)
   }  // for (int m = 0; m < numSub; m++)
   PrintParallelProfile("rom.assemble_local_domains",
                        MPI_Wtime() - profile_start, numSubLoc);

   // interface matrices.
   profile_start = MPI_Wtime();
   for (int p = 0; p < topol_handler->GetNumPorts(); p++)
   {
      const PortInfo *pInfo = topol_handler->GetPortInfo(p);
      const int p_type = topol_handler->GetPortType(p);
      MatrixBlocks *port_mat = rom_elems->port[p_type];

      const int m1 = pInfo->Mesh1;
      const int m2 = pInfo->Mesh2;

      int num_block = (separate_variable_basis) ? num_var : 1;

      Array<int> midx(0);

      for (int v = 0; v < num_block; v++)
         midx.Append(rom_handler->GetBlockIndex(m1, v));

      for (int v = 0; v < num_block; v++)
         midx.Append(rom_handler->GetBlockIndex(m2, v));

      AddToBlockMatrix(midx, midx, *port_mat, romMat);
   }
   PrintParallelProfile("rom.assemble_global_port_loop",
                        MPI_Wtime() - profile_start,
                        topol_handler->GetNumPorts());
}

void MultiBlockSolver::AssembleROMMat(SparseBlockMatrixBuilder &builder)
{
   assert(topol_mode == TopologyHandlerMode::COMPONENT);
   assert(rom_elems);

   const int ossub = topol_handler->GlobalSubdomainIndex(0);
   for (int mm = 0; mm < numSubLoc; mm++)
   {
      const int m = mm + ossub;
      const int m_global = topol_handler->GlobalSubdomainIndex(mm);
      const int c_type = topol_handler->GetMeshType(m_global);
      const int num_block = separate_variable_basis ? num_var : 1;

      Array<int> midx(num_block);
      for (int v = 0; v < num_block; v++)
         midx[v] = rom_handler->GetBlockIndex(m, v);

      AddToSparseBlockMatrix(midx, midx, *rom_elems->comp[c_type], builder);

      Array<int> *bdr_c2g = topol_handler->GetBdrAttrComponentToGlobalMap(mm);
      for (int b = 0; b < bdr_c2g->Size(); b++)
      {
         const int global_idx = global_bdr_attributes.Find((*bdr_c2g)[b]);
         if (global_idx < 0 || !BCExistsOnBdr(global_idx) ||
             bdr_type[global_idx] == BoundaryType::NEUMANN)
            continue;
         AddToSparseBlockMatrix(midx, midx,
                                *(*(rom_elems->bdr[c_type]))[b], builder);
      }
   }

   for (int p = 0; p < topol_handler->GetNumPorts(); p++)
   {
      const PortInfo *pInfo = topol_handler->GetPortInfo(p);
      const int p_type = topol_handler->GetPortType(p);
      const int num_block = separate_variable_basis ? num_var : 1;
      Array<int> midx(0);
      for (int v = 0; v < num_block; v++)
         midx.Append(rom_handler->GetBlockIndex(pInfo->Mesh1, v));
      for (int v = 0; v < num_block; v++)
         midx.Append(rom_handler->GetBlockIndex(pInfo->Mesh2, v));
      AddToSparseBlockMatrix(midx, midx, *rom_elems->port[p_type], builder);
   }
}

void MultiBlockSolver::InitVisualization(const std::string& output_path)
{
   if (!visual.save) return;

   std::string file_prefix;
   if (output_path != "")
      file_prefix = output_path;
   else
   {
      assert(visual.prefix != "");
      assert(visual.dir != "");
      file_prefix = visual.dir + "/" + visual.prefix;
   }

   if (visual.unified_view)
      InitUnifiedParaview(file_prefix);
   else
      InitIndividualParaview(file_prefix);
}

void MultiBlockSolver::InitIndividualParaview(const std::string& file_prefix)
{
   assert(var_names.size() == num_var);
   assert((visual.domain_offset >= 0) && (visual.domain_offset < numSubLoc));
   assert((visual.domain_interval > 0) && (visual.domain_interval <= numSubLoc));
   paraviewColls.SetSize(numSubLoc);
   paraviewColls = NULL;

   std::string error_type, tmp;
   if (visual.save_error)
   {
      error_visual.SetSize(num_var * numSubLoc);
      error_visual = NULL;
      for (int k = 0; k < error_visual.Size(); k++)
         error_visual[k] = new GridFunction(fes[k]);
      
      error_type = "_abs_error";
   }

   for (int m = 0; m < numSubLoc; m++) {
      if ((m < visual.domain_offset) || (m % visual.domain_interval != 0)) continue;

      ostringstream oss;
      // Each subdomain needs to be save separately.
      oss << file_prefix << "_" << std::to_string(m);

      paraviewColls[m] = new ParaViewDataCollection(oss.str().c_str(), &(*meshes[m]));
      paraviewColls[m]->SetLevelsOfDetail(order);
      paraviewColls[m]->SetHighOrderOutput(true);
      paraviewColls[m]->SetPrecision(8);

      for (int v = 0, idx = m * num_var; v < num_var; v++, idx++)
      {
         paraviewColls[m]->RegisterField(var_names[v].c_str(), us[idx]);
         if (visual.save_error)
         {
            tmp = var_names[v] + error_type;
            paraviewColls[m]->RegisterField(tmp.c_str(), error_visual[idx]);
         }
      }
      paraviewColls[m]->SetOwnData(false);
   }
}

void MultiBlockSolver::InitUnifiedParaview(const std::string& file_prefix)
{
   pmesh = topol_handler->GetGlobalMesh();

   global_fes.SetSize(num_var);
   global_us_visual.SetSize(num_var);
   for (int v = 0; v < num_var; v++)
   {
      global_fes[v] = new FiniteElementSpace(pmesh, fec[v], vdim[v]);
      global_us_visual[v] = new GridFunction(global_fes[v]);
   }

   std::string error_type, tmp;
   if (visual.save_error)
   {
      error_visual.SetSize(num_var * numSub);
      error_visual = NULL;
      for (int k = 0; k < error_visual.Size(); k++)
         error_visual[k] = new GridFunction(fes[k]);

      global_error_visual.SetSize(num_var);
      for (int v = 0; v < num_var; v++)
         global_error_visual[v] = new GridFunction(global_fes[v]);

      error_type = "_abs_error";
   }

   // TODO: For truly bottom-up case, when the parent mesh does not exist?
   mfem_warning("Paraview is unified. Any overlapped interface dof data will not be shown.\n");
   paraviewColls.SetSize(1);

   paraviewColls[0] = new ParaViewDataCollection(file_prefix.c_str(), pmesh);
   paraviewColls[0]->SetLevelsOfDetail(order);
   paraviewColls[0]->SetHighOrderOutput(true);
   paraviewColls[0]->SetPrecision(8);

   for (int v = 0; v < num_var; v++)
   {
      paraviewColls[0]->RegisterField(var_names[v].c_str(), global_us_visual[v]);
      if (visual.save_error)
      {
         tmp = var_names[v] + error_type;
         paraviewColls[0]->RegisterField(tmp.c_str(), global_error_visual[v]);
      }
   }
   paraviewColls[0]->SetOwnData(false);
}

void MultiBlockSolver::SaveVisualization()
{
   if (!visual.save) return;

   if (visual.unified_view)
   {
      mfem_warning("Paraview is unified. Any overlapped interface dof data will not be shown.\n");
      topol_handler->TransferToGlobal(us, global_us_visual, num_var);
      if (visual.save_error)
         topol_handler->TransferToGlobal(error_visual, global_error_visual, num_var);
   }
   else
   {
      assert((visual.domain_offset >= 0) && (visual.domain_offset < numSub));
      assert((visual.domain_interval > 0) && (visual.domain_interval <= numSub));
   }

   for (int m = 0; m < paraviewColls.Size(); m++)
   {
      if ((!visual.unified_view) && ((m < visual.domain_offset) || (m % visual.domain_interval != 0))) continue;
      assert(paraviewColls[m]);
      paraviewColls[m]->Save();
   }
};

void MultiBlockSolver::SaveVisualization(const int step, const double time)
{
   if (!visual.save) return;

   for (int m = 0; m < paraviewColls.Size(); m++)
   {
      paraviewColls[m]->SetCycle(step);
      paraviewColls[m]->SetTime(time);
   }

   SaveVisualization();
}

void MultiBlockSolver::SetParameterizedProblem(ParameterizedProblem *problem)
{
   assert(bdr_type.Size() == global_bdr_attributes.Size());
   for (int b = 0; b < global_bdr_attributes.Size(); b++)
   {
      int idx = problem->battr.Find(global_bdr_attributes[b]);
      if (idx < 0) continue;

      bdr_type[b] = problem->bdr_type[idx];
   }
}

void MultiBlockSolver::SaveSolution(std::string filename)
{
   if (!save_sol) return;

   if (filename == "")
   {
      filename = sol_dir + "/" + sol_prefix + ".h5";
   }
   printf("Saving the solution file %s ...", filename.c_str());

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   assert(file_id >= 0);

   // TODO: currently we only need solution vector. But we can add more data as we need.
   hdf5_utils::WriteDataset(file_id, "solution", *U);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
   printf("Done!\n");
}

void MultiBlockSolver::SaveSolutionWithTime(std::string filename, const int step, const double time)
{
   SaveSolution(filename);
   printf("Saving time/time step ...");

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   assert(file_id >= 0);

   // TODO: currently we only need solution vector. But we can add more data as we need.
   hdf5_utils::WriteAttribute(file_id, "timestep", step);
   hdf5_utils::WriteAttribute(file_id, "time", time);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
   printf("Done!\n");
}

void MultiBlockSolver::LoadSolution(const std::string &filename)
{
   // solution vector U must be instantiated beforehand.
   assert(U);
   printf("Loading the solution file %s ...", filename.c_str());
   if (!FileExists(filename))
      mfem_error("file does not exist!\n");
   

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);

   // Check the size of the solution in the file.
   hsize_t *dims = NULL;
   int ndims = hdf5_utils::GetDatasetSize(file_id, "solution", dims);
   bool size_match = (ndims == 1) && (dims[0] == U->Size());

   // Read the solution.
   if (size_match)
      hdf5_utils::ReadDataset(file_id, "solution", *U);

   // TODO: currently we only need solution vector. But we can add more data as we need.
   errf = H5Fclose(file_id);
   assert(errf >= 0);
   if (size_match)
      printf("Done!\n");
   else
      mfem_error("solution size does not match!\n");
}

void MultiBlockSolver::LoadSolutionWithTime(const std::string &filename, int &step, double &time)
{
   LoadSolution(filename);
   printf("Loading time/time step ...");

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);

   // TODO: currently we only need solution vector. But we can add more data as we need.
   hdf5_utils::ReadAttribute(file_id, "timestep", step);
   hdf5_utils::ReadAttribute(file_id, "time", time);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
   printf("Done!\n");
}

void MultiBlockSolver::CopySolution(BlockVector *input_sol)
{
   assert(input_sol->NumBlocks() == U->NumBlocks());
   for (int b = 0; b < U->NumBlocks(); b++)
      assert(input_sol->BlockSize(b) == U->BlockSize(b));
   *U = *input_sol;
}

void MultiBlockSolver::InitROMHandler()
{
   rom_handler = new MFEMROMHandler(topol_handler, var_offsets, var_names, separate_variable_basis);
   rom_handler->LoadBalanceROMBlocks(rank, nproc);

   /*
      Component-wise ROM assembly currently relies on the topology and ROM
      load balancers assigning the same contiguous subdomain range to each
      rank.  Make that implicit requirement explicit before the local-row
      matrix assembly path starts depending on it.

      A separate basis per variable creates multiple ROM blocks per
      subdomain, so that case needs a more general ownership mapping and is
      intentionally not covered by this check.
   */
   if ((topol_mode == TopologyHandlerMode::COMPONENT) &&
       !separate_variable_basis)
   {
      const std::array<int, 2> owned_blocks = rom_handler->GetLocalBlocks();
      MFEM_VERIFY(owned_blocks[1] - owned_blocks[0] == numSubLoc,
                  "ROM and topology load balancers assigned different local block counts");
      for (int m = 0; m < numSubLoc; m++)
      {
         MFEM_VERIFY(topol_handler->GlobalSubdomainIndex(m) == owned_blocks[0] + m,
                     "ROM and topology load balancers assigned different block ranges");
      }
   }

   if (!(topol_mode == TopologyHandlerMode::COMPONENT))
      return;

   GetComponentFESpaces(comp_fes);
   rom_elems = new ROMLinearElement(topol_handler, comp_fes, separate_variable_basis);
}

void MultiBlockSolver::GetBasisTags(std::vector<BasisTag> &basis_tags)
{
   if (separate_variable_basis)
   {
      basis_tags.resize(numSub * num_var);
      for (int m = 0, idx = 0; m < numSub; m++)
         for (int v = 0; v < num_var; v++, idx++)
            basis_tags[idx] = GetBasisTag(m, topol_handler, var_names[v]);
   }
   else
   {
      basis_tags.resize(numSub);
      for (int m = 0; m < numSub; m++)
         basis_tags[m] = GetBasisTag(m, topol_handler);
   }
}

BlockVector* MultiBlockSolver::PrepareSnapshots(std::vector<BasisTag> &basis_tags)
{
   BlockVector *U_snapshots = NULL;

   // View vector for U.
   if (separate_variable_basis)
      U_snapshots = new BlockVector(U->GetData(), var_offsets);
   else
      U_snapshots = new BlockVector(U->GetData(), domain_offsets);

   GetBasisTags(basis_tags);
   assert(U_snapshots->NumBlocks() == basis_tags.size());

   return U_snapshots;
}

void MultiBlockSolver::SaveSnapshots(SampleGenerator *sample_generator)
{
   assert(sample_generator);

   /* split the solution into each component with the corresponding tag */
   std::vector<BasisTag> basis_tags;
   BlockVector *U_snapshots = PrepareSnapshots(basis_tags);

   Array<int> col_idxs;
   sample_generator->SaveSnapshot(U_snapshots, basis_tags, col_idxs);
   sample_generator->SaveSnapshotPorts(topol_handler, col_idxs);

   /* delete only the view vector, not the data itself. */
   delete U_snapshots;
   return;
}

void MultiBlockSolver::ProjectRHSOnReducedBasis()
{
   // View vector for RHS.
   BlockVector *RHS_domain = NULL;
   if (separate_variable_basis)
      RHS_domain = new BlockVector(RHS->GetData(), var_offsets); 
   else
      RHS_domain = new BlockVector(RHS->GetData(), domain_offsets);
      
   rom_handler->ProjectRHSOnReducedBasis(RHS_domain);

   delete RHS_domain;
}

void MultiBlockSolver::SolveROM()
{
   // View vector for U.
   BlockVector *U_domain = NULL;
   if (separate_variable_basis)
      U_domain = new BlockVector(U->GetData(), var_offsets); 
   else
      U_domain = new BlockVector(U->GetData(), domain_offsets);
   
   rom_handler->Solve(U_domain);
}

void MultiBlockSolver::ComputeSubdomainErrorAndNorm(GridFunction *fom_sol, GridFunction *rom_sol, double &error, double &norm)
{
   assert(fom_sol && rom_sol);
   const int vec_dim = fom_sol->FESpace()->GetVDim();

   if (vec_dim == 1)
   {
      ConstantCoefficient zero(0.0);
      GridFunctionCoefficient rom_sol_coeff(rom_sol);

      norm = fom_sol->ComputeLpError(2, zero);
      error = fom_sol->ComputeLpError(2, rom_sol_coeff);
   }
   else
   {
      Vector zero_v(vec_dim);
      zero_v = 0.0;
      VectorConstantCoefficient zero(zero_v);
      VectorGridFunctionCoefficient rom_sol_coeff(rom_sol);

      norm = fom_sol->ComputeLpError(2, zero);
      error = fom_sol->ComputeLpError(2, rom_sol_coeff);
   }
}

void MultiBlockSolver::ComputeRelativeError(Array<GridFunction *> fom_sols, Array<GridFunction *> rom_sols, Vector &error)
{
   assert(fom_sols.Size() == (num_var * numSubLoc));
   assert(rom_sols.Size() == (num_var * numSubLoc));

   Vector norm(num_var);
   error.SetSize(num_var);
   norm = 0.0; error = 0.0;
   for (int m = 0, idx = 0; m < numSubLoc; m++)
   {
      for (int v = 0; v < num_var; v++, idx++)
      {
         assert(fom_sols[idx] && rom_sols[idx]);
         double var_norm = 0.0, var_error = 0.0;
         ComputeSubdomainErrorAndNorm(fom_sols[idx], rom_sols[idx], var_error, var_norm);
         norm[v] += var_norm * var_norm;
         error[v] += var_error * var_error;
      }
   }

   for (int v = 0; v < num_var; v++)
   {
      double sum = 0.0;
      MPI_Allreduce(&norm[v], &sum, 1, MPI_DOUBLE, MPI_SUM, topol_handler->GetComm());
      norm[v] = sqrt(sum);

      sum = 0.0;
      MPI_Allreduce(&error[v], &sum, 1, MPI_DOUBLE, MPI_SUM, topol_handler->GetComm());
      error[v] = sqrt(sum);

      if (rank == 0) printf("Variable %d global absolute error: %.5E\n", v, error[v]);
      error[v] /= norm[v];
      if (rank == 0) printf("Variable %d global relative error: %.5E\n", v, error[v]);
   }
}

void MultiBlockSolver::CompareSolution(BlockVector &test_U, Vector &error)
{
   assert(test_U.NumBlocks() == U->NumBlocks());
   for (int b = 0; b < U->NumBlocks(); b++)
      assert(test_U.BlockSize(b) == U->BlockSize(b));

   Array<GridFunction *> test_us;
   test_us.SetSize(num_var * numSubLoc);
   for (int k = 0; k < test_us.Size(); k++)
   {
      test_us[k] = new GridFunction(fes[k], test_U.GetBlock(k), 0);

      // BC's are weakly constrained and there is no essential dofs.
      // Does this make any difference?
      test_us[k]->SetTrueVector();
   }

   // Compare the solution.
   // Maximum L2-error / L2-norm over all variables.
   error.SetSize(num_var);
   ComputeRelativeError(us, test_us, error);

   if (visual.save_error)
      for (int k = 0; k < test_us.Size(); k++)
         subtract(*test_us[k], *us[k], *error_visual[k]);

   DeletePointers(test_us);
}
