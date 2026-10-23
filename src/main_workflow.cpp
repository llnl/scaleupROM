// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "main_workflow.hpp"
#include "component_topology_handler.hpp"
#include "multiblock_solver.hpp"
#include "poisson_solver.hpp"
#include "linelast_solver.hpp"
#include "stokes_solver.hpp"
#include "steady_ns_solver.hpp"
#include "advdiff_solver.hpp"
#include "unsteady_ns_solver.hpp"
#include "etc.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

double dbc2(const Vector &x, double t)
{
   return 0.1 - 0.1 * (x(1) - 1.0) * (x(1) - 1.0);
}

double dbc4(const Vector &x, double t)
{
   return -0.1 + 0.1 * (x(1) - 1.0) * (x(1) - 1.0);
}

void RunExample()
{
   PoissonSolver test;

   test.InitVariables();
   test.InitVisualization();

   test.AddBCFunction(dbc2, 2);
   test.AddBCFunction(dbc4, 4);
   test.AddRHSFunction(1.0);

   test.BuildOperators();

   test.SetupBCOperators();

   test.Assemble();

   test.Solve();
   test.SaveVisualization();
}

MultiBlockSolver* InitSolver()
{
   std::string solver_type = config.GetRequiredOption<std::string>("main/solver");
   MultiBlockSolver *solver = NULL;
   if (solver_type == "poisson")       { solver = new PoissonSolver; }
   else if (solver_type == "stokes")   { solver = new StokesSolver; }
   else if (solver_type == "steady-ns")   { solver = new SteadyNSSolver; }
   else if (solver_type == "linelast")   { solver = new LinElastSolver; }
   else if (solver_type == "adv-diff")   { solver = new AdvDiffSolver; }
   else if (solver_type == "unsteady-ns")   { solver = new UnsteadyNSSolver; }
   else
   {
      printf("Unknown MultiBlockSolver %s!\n", solver_type.c_str());
      exit(-1);
   }

   return solver;
}

SampleGenerator* InitSampleGenerator(MPI_Comm comm)
{
   SampleGenerator* generator = NULL;

   std::string type = config.GetOption<std::string>("sample_generation/type", "base");

   if (type == "base")
   {
      generator = new SampleGenerator(comm);
   }
   else if (type == "random")
   {
      generator = new RandomSampleGenerator(comm);
   }
   else
   {
      mfem_error("Unknown sample generator type!\n");
   }

   return generator;
}

std::vector<BasisTag> GetGlobalBasisTagList(const TopologyHandlerMode &topol_mode, bool separate_variable_basis)
{
   std::vector<BasisTag> basis_tags(0);

   std::vector<std::string> component_list(0);
   if (topol_mode == TopologyHandlerMode::SUBMESH)
   {
      TopologyHandler *topol_handler = new SubMeshTopologyHandler();
      for (int c = 0; c < topol_handler->GetNumComponents(); c++)
         component_list.push_back(topol_handler->GetComponentName(c));

      delete topol_handler;
   }
   else if (topol_mode == TopologyHandlerMode::COMPONENT)
   {
      YAML::Node component_dict = config.FindNode("mesh/component-wise/components");
      assert(component_dict);
      for (int p = 0; p < component_dict.size(); p++)
      component_list.push_back(config.GetRequiredOptionFromDict<std::string>("name", component_dict[p]));
   }
   else
      mfem_error("GetGlobalBasisTagList - TopologyHandlerMode is not set!\n");

   std::vector<std::string> var_list(0);
   if (separate_variable_basis)
   {
      std::string solver_type = config.GetRequiredOption<std::string>("main/solver");
      if (solver_type == "poisson")          var_list = PoissonSolver::GetVariableNames();
      else if (solver_type == "stokes")      var_list = StokesSolver::GetVariableNames();
      else if (solver_type == "steady-ns")   var_list = SteadyNSSolver::GetVariableNames();
      else if (solver_type == "linelast")   var_list = LinElastSolver::GetVariableNames();
      else if (solver_type == "unsteady-ns")   var_list = UnsteadyNSSolver::GetVariableNames();
      else
      {
         printf("Unknown MultiBlockSolver %s!\n", solver_type.c_str());
         exit(-1);
      }
   }

   for (int c = 0; c < component_list.size(); c++)
   {
      if (separate_variable_basis)
         for (int v = 0; v < var_list.size(); v++)
            basis_tags.push_back(BasisTag(component_list[c], var_list[v]));
      else
         basis_tags.push_back(BasisTag(component_list[c]));
   }

   return basis_tags;
}

void GenerateSamples(MPI_Comm comm)
{
   // save the original config.dict_
   YAML::Node dict0 = YAML::Clone(config.dict_);
   ParameterizedProblem *problem = InitParameterizedProblem();
   SampleGenerator *sample_generator = InitSampleGenerator(comm);
   SampleGeneratorType sample_gen_type = sample_generator->GetType();
   sample_generator->SetParamSpaceSizes();
   MultiBlockSolver *test = NULL;

   int s = 0;
   while (s < sample_generator->GetTotalSampleSize())
   {
      if (!sample_generator->IsMyJob(s)) {
	s++;
	continue;
      }

      // NOTE: this will change config.dict_
      sample_generator->SetSampleParams(s);

      test = InitSolver();
      test->InitVariables();
      if (test->UseRom())
         test->InitROMHandler();

      problem->SetSingleRun();
      test->SetParameterizedProblem(problem);

      int file_idx = s + sample_generator->GetFileOffset();
      const std::string visual_path = sample_generator->GetSamplePath(file_idx, test->GetVisualizationPrefix());
      std::string sol_file = sample_generator->GetSamplePath(file_idx, test->GetSolutionFilePrefix());
      sol_file += ".h5";
      test->InitVisualization(visual_path);
      test->BuildOperators();
      test->SetupBCOperators();
      test->Assemble();

      bool converged = test->Solve(sample_generator);
      if (!converged)
      {
         // If deterministic, terminate the sampling here.
         if (sample_gen_type == BASE)
            mfem_error("A sample solution fails to converge!\n");
         else if (sample_gen_type == RANDOM)
         {
            // if random, try another sample.
            mfem_warning("A sample solution failed to converge. Trying another sample.\n");
            delete test;
            continue;
         }
      }
         
      test->SaveSolution(sol_file);
      test->SaveVisualization();

      sample_generator->ReportStatus(s);

      delete test;

      s++;
   }
   sample_generator->WriteSnapshots();
   sample_generator->WriteSnapshotPorts();

   delete sample_generator;
   delete problem;
   // restore the original config.dict_
   config.dict_ = dict0;
}

void CollectSamples(SampleGenerator *sample_generator)
{
   std::string mode = config.GetOption<std::string>("sample_collection/mode", "basis");
   std::string basis_prefix = config.GetOption<std::string>("basis/prefix", "basis");

   if (mode == "basis")
      CollectSamplesByBasis(sample_generator, basis_prefix);
   else if (mode == "port")
      CollectSamplesByPort(sample_generator, basis_prefix);
   else
      mfem_error("CollectSamples: unknown sample collection mode!\n");
}

void CollectSamplesByPort(SampleGenerator *sample_generator, const std::string &basis_prefix)
{
   // parse the sample snapshot file list.
   std::vector<std::string> file_list = config.GetOption<std::vector<std::string>>(
                                 "sample_collection/port_files", std::vector<std::string>(0));
   YAML::Node port_format = config.FindNode("sample_collection/port_fileformat");
   // if file list is specified with a format, parse through the format.
   if (port_format)
   {
      FilenameParam port_param("", port_format);
      port_param.ParseFilenames(file_list);
   }

   // if additional inputs are not specified for port files, set default port file name.
   if (file_list.size() == 0)
      file_list.push_back(sample_generator->GetSamplePrefix() + ".port.h5");

   for (int f = 0; f < file_list.size(); f++)
      sample_generator->CollectSnapshotsByPort(basis_prefix, file_list[f]);
}

void CollectSamplesByBasis(SampleGenerator *sample_generator, const std::string &basis_prefix)
{
   assert(sample_generator);

   TopologyHandlerMode topol_mode = SetTopologyHandlerMode();
   bool separate_variable_basis = config.GetOption<bool>("model_reduction/separate_variable_basis", false);

   // Find the all required basis tags.
   std::vector<BasisTag> basis_tags = GetGlobalBasisTagList(topol_mode, separate_variable_basis);

   // tag-specific optional inputs.
   YAML::Node basis_list = config.FindNode("basis/tags");

   // loop over the required basis tag list.
   for (int p = 0; p < basis_tags.size(); p++)
   {
      std::vector<std::string> file_list(0);
      std::string default_filename = sample_generator->GetBaseFilename(sample_generator->GetSamplePrefix(), basis_tags[p]);
      default_filename += "_snapshot";
      FindSnapshotFilesForBasis(basis_tags[p], default_filename, file_list);
      assert(file_list.size() > 0);

      sample_generator->CollectSnapshotsByBasis(basis_prefix, basis_tags[p], file_list);
   }  // for (int p = 0; p < basis_tags.size(); p++)
}

void TrainROM(MPI_Comm comm)
{
   SampleGenerator *sample_generator = InitSampleGenerator(comm);

   std::string basis_prefix = config.GetOption<std::string>("basis/prefix", "basis");
   CollectSamples(sample_generator);

   sample_generator->FormReducedBasis(basis_prefix);

   AuxiliaryTrainROM(comm);

   delete sample_generator;
}

void AuxiliaryTrainROM(MPI_Comm comm)
{
   std::string solver_type = config.GetRequiredOption<std::string>("main/solver");
   bool separate_variable_basis = config.GetOption<bool>("model_reduction/separate_variable_basis", false);

   /* Supremizer enrichment */
   if ((separate_variable_basis) &&
       ((solver_type == "stokes") || (solver_type == "steady-ns") || (solver_type == "unsteady-ns")))
   {
      ParameterizedProblem *problem = InitParameterizedProblem();
      StokesSolver *solver = NULL;
      if (solver_type == "stokes")           { solver = new StokesSolver; }
      else if (solver_type == "steady-ns")   { solver = new SteadyNSSolver; }
      else if (solver_type == "unsteady-ns")   { solver = new UnsteadyNSSolver; }

      if (!solver->UseRom()) mfem_error("ROM must be enabled for supremizer enrichment!\n");

      solver->InitVariables();
      solver->InitROMHandler();
      // This time needs to be ROMHandler, in order not to run StokesSolver::LoadSupremizer.
      solver->GetROMHandler()->LoadReducedBasis();

      solver->EnrichSupremizer();

      delete problem;
      delete solver;
   }
}

void TrainEQP(MPI_Comm comm)
{
   SampleGenerator *sample_generator = InitSampleGenerator(comm);

   std::string basis_prefix = config.GetOption<std::string>("basis/prefix", "basis");
   CollectSamples(sample_generator);

   /* EQP NNLS procedure */
   std::string eqp_str = config.GetOption<std::string>("model_reduction/nonlinear_handling", "none");
   if (eqp_str != "eqp")
      mfem_error("ROM nonlinear handling is not eqp!\n");

   MultiBlockSolver *test = NULL;
   test = InitSolver();
   test->InitVariables();
   test->InitROMHandler();

   if (!test->IsNonlinear())
   {
      mfem_warning("The physics solver is not nonlinear. Exiting TrainEQP.\n");
      delete test;
      return;
   }

   if (!test->UseRom()) mfem_error("ROM must be enabled for EQP training!\n");

   test->LoadReducedBasis();
   test->AllocateROMNlinElems();

   ROMHandlerBase *rom = test->GetROMHandler();
   ROMBuildingLevel save_operator = rom->GetBuildingLevel();
   TopologyHandlerMode topol_mode = test->GetTopologyMode();

   if (topol_mode == TopologyHandlerMode::SUBMESH)
      printf("using SubMesh topology.\n");
   else if (topol_mode == TopologyHandlerMode::COMPONENT)
      printf("using Component-wise topology.\n");
   else
      mfem_error("Unknown TopologyHandler Mode!\n");

   std::string oper_prefix = rom->GetOperatorPrefix();
   switch (save_operator)
   {
      case ROMBuildingLevel::COMPONENT:
      {
         if (topol_mode == TopologyHandlerMode::SUBMESH)
            mfem_error("Submesh does not support component rom building level!\n");

         test->TrainROMEQPElems(sample_generator);
         test->SaveROMNlinElems(oper_prefix);
         break;
      }
      case ROMBuildingLevel::GLOBAL:
      {
         mfem_error("TrainEQP: not implemented for global yet!\n");
         // test->TrainEQP(sample_generator);
         // test->SaveEQP();
         break;
      }
      case ROMBuildingLevel::NONE:
      default:
         mfem_error("TrainEQP: save_operator level must be either component or global!\n");
         break;
   }

   delete test;
}

void FindSnapshotFilesForBasis(const BasisTag &basis_tag, const std::string &default_filename, std::vector<std::string> &file_list)
{
   file_list.clear();

   // tag-specific optional inputs.
   YAML::Node basis_list = config.FindNode("basis/tags");

   // if optional inputs are specified, parse them first.
   if (basis_list)
   {
      // Find if additional inputs are specified for basis_tag.
      YAML::Node basis_tag_input = config.LookUpFromDict("name", basis_tag.print(), basis_list);
      
      // If basis_tag has additional inputs, parse them.
      if (basis_tag_input)
      {
         // parse the sample snapshot file list.
         file_list = config.GetOptionFromDict<std::vector<std::string>>(
                     "snapshot_files", std::vector<std::string>(0), basis_tag_input);
         YAML::Node snapshot_format = config.FindNodeFromDict("snapshot_format", basis_tag_input);
         // if file list is specified with a format, parse through the format.
         if (snapshot_format)
         {
            FilenameParam snapshot_param("", snapshot_format);
            snapshot_param.ParseFilenames(file_list);
         }
      }  // if (basis_tag_input)
   }

   // if additional inputs are not specified for snapshot files, set default snapshot file name.
   if (file_list.size() == 0)
      file_list.push_back(default_filename);
}

void BuildROM(MPI_Comm comm)
{
   ParameterizedProblem *problem = InitParameterizedProblem();
   MultiBlockSolver *test = NULL;

   test = InitSolver();
   if (!test->UseRom()) mfem_error("ROM must be enabled for BuildROM!\n");
   test->InitVariables();
   test->InitROMHandler();
   // test->InitVisualization();

   // The ROM operator will be built based on the parameter specified for single-run.
   problem->SetSingleRun();
   test->SetParameterizedProblem(problem);

   // TODO: there are skippable operations depending on rom/fom mode.
   test->BuildOperators();
   test->SetupBCOperators();
   test->LoadReducedBasis();

   if (test->IsNonlinear())
      test->AllocateROMNlinElems();
   
   TopologyHandlerMode topol_mode = test->GetTopologyMode();
   ROMHandlerBase *rom = test->GetROMHandler();
   ROMBuildingLevel save_operator = rom->GetBuildingLevel();

   // NOTE(kevin): global operator required only for global rom operator.
   if (save_operator == ROMBuildingLevel::GLOBAL)
      test->Assemble();

   std::string oper_prefix = rom->GetOperatorPrefix();
   switch (save_operator)
   {
      case ROMBuildingLevel::COMPONENT:
      {
         if (topol_mode == TopologyHandlerMode::SUBMESH)
            mfem_error("Submesh does not support component rom building level!\n");

         test->BuildROMLinElems();
         test->SaveROMLinElems(oper_prefix + ".h5");

         if ((test->IsNonlinear()) && (rom->GetNonlinearHandling() == NonlinearHandling::TENSOR))
         {
            test->BuildROMTensorElems();
            test->SaveROMNlinElems(oper_prefix);
         }
         break;
      }
      case ROMBuildingLevel::GLOBAL:
      {
         test->ProjectOperatorOnReducedBasis();
         test->SaveROMOperator(oper_prefix + ".h5");
         break;
      }
      case ROMBuildingLevel::NONE:
      {
         printf("BuildROM - ROM building level is set to none. No ROM is saved.\n");
         break;
      }
   }  // switch (save_operator)

   test->SaveBasisVisualization();

   delete test;
   delete problem;
}

double SingleRun(MPI_Comm comm, const std::string output_file)
{
   double profile_start = MPI_Wtime();
   if (config.GetOption<bool>("single_run/choose_from_random_sample", false))
   {
      RandomSampleGenerator *generator = new RandomSampleGenerator(comm);
      generator->SetParamSpaceSizes();
      int idx = UniformRandom(0, generator->GetTotalSampleSize()-1);
      // NOTE: this will change config.dict_
      generator->SetSampleParams(idx);
      delete generator;
   }

   ParameterizedProblem *problem = InitParameterizedProblem();
   PrintParallelProfile("init.parameterized_problem",
                        MPI_Wtime() - profile_start);
   profile_start = MPI_Wtime();
   MultiBlockSolver *test = InitSolver();
   PrintParallelProfile("init.solver_and_topology", MPI_Wtime() - profile_start);
   profile_start = MPI_Wtime();
   test->InitVariables();
   if (test->UseRom()) test->InitROMHandler();
   test->InitVisualization();
   PrintParallelProfile("init.variables_rom_visualization",
                        MPI_Wtime() - profile_start);

   StopWatch solveTimer;
   std::string solveType = (test->UseRom()) ? "ROM" : "FOM";

   problem->SetSingleRun();
   test->SetParameterizedProblem(problem);

   // TODO: there are skippable operations depending on rom/fom mode.
   profile_start = MPI_Wtime();
   test->BuildRHSOperators();
   test->SetupRHSBCOperators();
   test->AssembleRHS();
   PrintParallelProfile("setup.rhs_operators_bc_assembly",
                        MPI_Wtime() - profile_start);

   const int num_var = test->GetNumVar();
   Vector rom_assemble(1), rom_solve(1), fom_assemble(1), fom_solve(1), error(num_var);
   rom_assemble = -1.0; rom_solve = -1.0;
   fom_assemble = -1.0; fom_solve = -1.0;
   error = -1.0;

   ROMHandlerBase *rom = NULL;
   if (test->UseRom())
   {
      rom = test->GetROMHandler();
      profile_start = MPI_Wtime();
      test->LoadReducedBasis();
      PrintParallelProfile("rom.load_reduced_basis",
                           MPI_Wtime() - profile_start);

      if (test->IsNonlinear())
         test->AllocateROMNlinElems();
   }

   solveTimer.Start();
   if (test->UseRom())
   {
      printf("ROM with ");
      ROMBuildingLevel save_operator = rom->GetBuildingLevel();
      TopologyHandlerMode topol_mode = test->GetTopologyMode();

      if (topol_mode == TopologyHandlerMode::SUBMESH)
         printf("using SubMesh topology.\n");
      else if (topol_mode == TopologyHandlerMode::COMPONENT)
         printf("using Component-wise topology.\n");
      else
         mfem_error("Unknown TopologyHandler Mode!\n");

      std::string filename = rom->GetOperatorPrefix() + ".h5";
      if (save_operator == ROMBuildingLevel::COMPONENT)
      {
         if (topol_mode == TopologyHandlerMode::SUBMESH)
            mfem_error("Submesh does not support component rom building level!\n");

         printf("Loading ROM projected elements.. ");
         profile_start = MPI_Wtime();
         test->LoadROMLinElems(filename);
         PrintParallelProfile("rom.load_component_operators",
                              MPI_Wtime() - profile_start);
         printf("Done!\n");

         printf("Assembling ROM linear matrix.. ");
         profile_start = MPI_Wtime();
         test->AssembleROMMat();
         PrintParallelProfile("rom.assemble_total",
                              MPI_Wtime() - profile_start);
         printf("Done!\n");

         if (test->IsNonlinear())
         {
            test->LoadROMNlinElems(rom->GetOperatorPrefix());
            test->AssembleROMNlinOper();
         }
      }  // if (save_operator == ROMBuildingLevel::COMPONENT)
      else if (save_operator == ROMBuildingLevel::GLOBAL)
      {
         printf("Loading global operator file.. ");
         test->LoadROMOperatorFromFile(filename);
         printf("Done!\n");
      }  // if (save_operator == ROMBuildingLevel::GLOBAL)
      else if (save_operator == ROMBuildingLevel::NONE)
      {
         printf("Building operator file all the way from FOM.. ");
         test->BuildDomainOperators();
         test->SetupDomainBCOperators();
         test->AssembleOperator();
         test->ProjectOperatorOnReducedBasis();
         printf("Done!\n");
      }  // if (save_operator == ROMBuildingLevel::NONE)
      else
         mfem_error("SingleRun - Unknown ROMBuildingLevel!\n");

      printf("Projecting RHS to ROM.. ");
      profile_start = MPI_Wtime();
      test->ProjectRHSOnReducedBasis();
      PrintParallelProfile("rom.project_rhs", MPI_Wtime() - profile_start);
      printf("Done!\n");
   }  // if (test->UseRom())
   else
   {
      test->BuildDomainOperators();
      test->SetupDomainBCOperators();
      test->AssembleOperator();
   }  // not if (test->UseRom())
   solveTimer.Stop();
   printf("%s-assemble time: %f seconds.\n", solveType.c_str(), solveTimer.RealTime());

   if (test->UseRom())
      rom_assemble = solveTimer.RealTime();
   else
      fom_assemble = solveTimer.RealTime();

   solveTimer.Clear();
   solveTimer.Start();
   if (test->UseRom())
   {
      test->SolveROM();
   }
   else
   {
      // TODO: move matrix assembly to here.
      test->Solve();
   }
   solveTimer.Stop();
   printf("%s-solve time: %f seconds.\n", solveType.c_str(), solveTimer.RealTime());

   if (test->UseRom())
      rom_solve = solveTimer.RealTime();
   else
      fom_solve = solveTimer.RealTime();

   profile_start = MPI_Wtime();
   /* save the ROM system for analysis/debug */
   bool save_rom = config.GetOption<bool>("model_reduction/save_linear_system/enabled", false);
   if (save_rom)
   {
      std::string rom_prefix = config.GetRequiredOption<std::string>("model_reduction/save_linear_system/prefix");
      rom->SaveRomSystem(rom_prefix);
   }

   bool compare_sol = config.GetOption<bool>("model_reduction/compare_solution/enabled", false);
   bool load_sol = config.GetOption<bool>("model_reduction/compare_solution/load_solution", false);
   if (test->UseRom() && compare_sol)
   {
      BlockVector *romU = test->GetSolutionCopy();

      if (load_sol)
      {
         printf("Comparing with the existing FOM solution.\n");
         std::string fom_file = config.GetRequiredOption<std::string>("model_reduction/compare_solution/fom_solution_file");
         test->LoadSolution(fom_file);
      }
      else
      {
         solveTimer.Clear();
         solveTimer.Start();
         test->BuildDomainOperators();
         test->SetupDomainBCOperators();
         test->AssembleOperator();
         solveTimer.Stop();
         printf("FOM-assembly time: %f seconds.\n", solveTimer.RealTime());
         fom_assemble = solveTimer.RealTime();

         solveTimer.Clear();
         solveTimer.Start();
         test->Solve();
         solveTimer.Stop();
         printf("FOM-solve time: %f seconds.\n", solveTimer.RealTime());
         fom_solve = solveTimer.RealTime();
      }

      test->CompareSolution(*romU, error);

      bool save_reduced_sol = config.GetOption<bool>("model_reduction/compare_solution/save_reduced_solution", false);
      if (save_reduced_sol)
      {
         ROMHandlerBase *rom = test->GetROMHandler();
         rom->SaveReducedSolution("rom_reduced_sol.txt");

         // use ROMHandler::reduced_rhs as a temporary variable.
         rom->ProjectRHSOnReducedBasis(test->GetSolution());
         rom->SaveReducedRHS("fom_reduced_sol.txt");
      }

      // Recover the original ROM solution.
      test->CopySolution(romU);

      delete romU;
   }

   // save results to output file.
   if (output_file.length() > 0)
   {
      hid_t file_id;
      herr_t errf = 0;
      file_id = H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      assert(file_id >= 0);

      hdf5_utils::WriteDataset(file_id, "rom_assemble", rom_assemble);
      hdf5_utils::WriteDataset(file_id, "rom_solve", rom_solve);
      hdf5_utils::WriteDataset(file_id, "fom_assemble", fom_assemble);
      hdf5_utils::WriteDataset(file_id, "fom_solve", fom_solve);
      hdf5_utils::WriteDataset(file_id, "rel_error", error);

      errf = H5Fclose(file_id);
      assert(errf >= 0);
   }

   // Save solution and visualization.
   test->SaveSolution();
   test->SaveVisualization();
   PrintParallelProfile("postprocess.output_and_comparison",
                        MPI_Wtime() - profile_start);
   
   profile_start = MPI_Wtime();
   delete test;
   delete problem;
   PrintParallelProfile("cleanup.solver_and_problem",
                        MPI_Wtime() - profile_start);

   // return the maximum error over all variables.
   return error.Max();
}

void PrintEQPCoords(MPI_Comm comm)
{
   ParameterizedProblem *problem = InitParameterizedProblem();
   MultiBlockSolver *test = InitSolver();
   test->InitVariables();

   assert(test->UseRom());
   if (!test->IsNonlinear())
   {
      printf("PrintEQPPoints: given solver is not nonlinear. Exiting\n");
      return;
   }

   test->InitROMHandler();

   StopWatch solveTimer;

   problem->SetSingleRun();
   test->SetParameterizedProblem(problem);

   // TODO: there are skippable operations depending on rom/fom mode.
   test->BuildRHSOperators();
   test->SetupRHSBCOperators();
   test->AssembleRHS();

   const int num_var = test->GetNumVar();

   ROMHandlerBase *rom = test->GetROMHandler();
   test->LoadReducedBasis();
   test->AllocateROMNlinElems();

   ROMBuildingLevel save_operator = rom->GetBuildingLevel();
   TopologyHandlerMode topol_mode = test->GetTopologyMode();
   assert(save_operator == ROMBuildingLevel::COMPONENT);
   assert(topol_mode != TopologyHandlerMode::SUBMESH);
   assert(rom->GetNonlinearHandling() == NonlinearHandling::EQP);

   test->LoadROMNlinElems(rom->GetOperatorPrefix());

   const std::string filename = config.GetRequiredOption<std::string>("model_reduction/eqp/save_coords");
   test->SaveEQPCoords(filename);
   
   delete test;
   delete problem;

   // return the maximum error over all variables.
   return;
}
