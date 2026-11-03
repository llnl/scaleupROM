// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef SCALEUPROM_MULTIBLOCK_SOLVER_HPP
#define SCALEUPROM_MULTIBLOCK_SOLVER_HPP

#include "topology_handler.hpp"
#include "interfaceinteg.hpp"
#include "interface_form.hpp"
#include "mfem.hpp"
#include "parameterized_problem.hpp"
#include "rom_handler.hpp"
#include "hdf5_utils.hpp"
#include "sample_generator.hpp"
#include "rom_element_collection.hpp"

// By convention we only use mfem namespace as default, not CAROM.
using namespace mfem;

class MultiBlockSolver
{

friend class ParameterizedProblem;

protected:
   /* MPI variables */
   int nproc;
   int rank;

   /*
      Base variables needed for all systems (potentially)
   */
   int order = 1;

   bool full_dg = true;
   int skip_zeros = 1;
   bool nonlinear_mode = false;

   TopologyHandlerMode topol_mode = NUM_TOPOL_MODE;
   TopologyHandler *topol_handler = NULL;

   // MultiBlockSolver does not own these. Owned by TopologyHandler.
   Mesh *pmesh = NULL;  // parent mesh. only available from SubMeshTopologyHandler.
   Array<Mesh*> meshes;

   // Information received from Topology Handler.
   int numSub;   // number of subdomains.
   int numSubLoc;   // number of local subdomains.
   int numSubStored;   // number of subdomains stored on this MPI rank.
   int dim;      // Spatial dimension.
   Array<int> global_bdr_attributes;   // boundary attributes of global system.

   // MPI partitioning data
   Array<int> localNumBlocks;
   std::array<int,2> localBlocks;
   Array<int> neighbors;

   // Solution dimension, by default -1 (scalar).
   int udim = -1;       // vector dimension of the entire solution variable
   int num_var = -1;    // number of variables
   Array<int> vdim;     // vector dimension of each variable
   std::vector<std::string> var_names;

   Array<int> block_offsets;  // Size(numSub * udim + 1). each block corresponds to a component of vector solution.
   Array<int> var_offsets; // Size(num_var * numSub + 1). each block corresponds to the variable of the solution in each domain.
   Array<int> domain_offsets; // Size(numSub + 1). each block corresponds to the vector solution in each domain.
   Array<int> num_vdofs;       // Size(numSub). number of vdofs of the vector solution in each subdomain.
   BlockVector *U, *RHS;

   // Each Variable is separated by distinct FiniteElementSpace.
   Array<FiniteElementCollection *> fec;  // Size(num_var);
   Array<FiniteElementSpace *> fes;       // Size(num_var * numSub);
   Array<GridFunction *> us;              // Size(num_var * numSub);

   int max_bdr_attr;
   int numBdr;
   Array<Array<int> *> bdr_markers;
   Array<BoundaryType> bdr_type; // Boundary condition types of (numBdr) array

   // MFEM solver options
   bool use_amg;
   bool direct_solve = false;

   // Saving solution in single run
   bool save_sol = false;
   std::string sol_dir = ".";
   std::string sol_prefix;

   /* visualization options */
   struct VisualizationOption {
      bool save = false;
      bool unified_view = false;
      std::string dir = ".";
      std::string prefix;
      /* for individual_paraview, support partial visualization. */
      int domain_offset = -1;
      int domain_interval = -1;

      /* time-dependent options */
      int time_interval = -1;

      /* visualizing rom solution error */
      bool save_error = false;
   } visual;

   // visualization variables
   Array<ParaViewDataCollection *> paraviewColls;
   // Used only for the unified visualization. Size(num_var).
   Array<FiniteElementSpace *> global_fes;
   Array<GridFunction *> global_us_visual;

   Array<GridFunction *> error_visual;    // Size(num_var * numSub)
   Array<GridFunction *> global_error_visual;   // point-wise error visualization

   // rom variables.
   ROMHandlerBase *rom_handler = NULL;
   bool use_rom = false;
   bool separate_variable_basis = false;

   // Used for bottom-up building, only with ComponentTopologyHandler.
   Array<FiniteElementSpace *> comp_fes;
   ROMLinearElement *rom_elems = NULL;

public:
   MultiBlockSolver();

   virtual ~MultiBlockSolver();

   // Parse some base input options. 
   void ParseInputs();

   // access
   const int GetDim() const { return dim; }
   const int GetNumSubdomains() const { return numSub; }
   const int GetNumBdr() const { return numBdr; }
   const int GetNumVar() const { return num_var; }
   Mesh* GetMesh(const int k) { return &(*meshes[k]); }
   GridFunction* GetGridFunction(const int k) { return us[k]; }
   const int GetDiscretizationOrder() const { return order; }
   const bool IsNonlinear() const { return nonlinear_mode; }
   const bool UseRom() const { return use_rom; }
   ROMHandlerBase* GetROMHandler() const { return rom_handler; }
   TopologyHandler* GetTopologyHandler() const { return topol_handler; }
   const bool IsVisualizationSaved() const { return visual.save; }
   const std::string GetSolutionFilePrefix() const { return sol_prefix; }
   const std::string GetVisualizationPrefix() const { return visual.prefix; }
   const TopologyHandlerMode GetTopologyMode() const { return topol_mode; }
   ParaViewDataCollection* GetParaViewColl(const int &k) { return paraviewColls[k]; }
   BlockVector* GetSolution() { return U; }
   BlockVector* GetSolutionCopy() { return new BlockVector(*U); }

   void SetSolutionSaveMode(const bool save_sol_);

   void GetVariableVector(const int &var_idx, BlockVector &global, BlockVector &var);
   void SetVariableVector(const int &var_idx, BlockVector &var, BlockVector &global);

   void SortBySubdomains(BlockVector &by_var, BlockVector &by_sub);
   void SortByVariables(BlockVector &by_sub, BlockVector &by_var);

   virtual void SetupBCVariables();
   virtual void AddBCFunction(std::function<double(const Vector &, double)> F, const int battr = -1)
   { mfem_error("Abstract method MultiBlockSolver::AddBCFunction!\n"); }
   virtual void AddBCFunction(const double &F, const int battr = -1)
   { mfem_error("Abstract method MultiBlockSolver::AddBCFunction!\n"); }
   virtual void AddBCFunction(std::function<void(const Vector &, double, Vector &)> F, const int battr = -1)
   { mfem_error("Abstract method MultiBlockSolver::AddBCFunction!\n"); }
   virtual void AddBCFunction(const Vector &F, const int battr = -1)
   { mfem_error("Abstract method MultiBlockSolver::AddBCFunction!\n"); }
   virtual void InitVariables() = 0;

   virtual void BuildOperators() = 0;
   virtual void BuildRHSOperators() = 0;
   virtual void BuildDomainOperators() = 0;
   
   void SetBdrType(const BoundaryType type, const int &global_battr_idx=-1);
   virtual bool BCExistsOnBdr(const int &global_battr_idx) = 0;
   virtual void SetupBCOperators() = 0;
   virtual void SetupRHSBCOperators() = 0;
   virtual void SetupDomainBCOperators() = 0;

   virtual void AddRHSFunction(std::function<double(const Vector &)> F)
   { mfem_error("Abstract method MultiBlockSolver::AddRHSFunction!\n"); }
   virtual void AddRHSFunction(const double F)
   { mfem_error("Abstract method MultiBlockSolver::AddRHSFunction!\n"); }
   virtual void AddRHSFunction(std::function<void(const Vector &, Vector &)> F)
   { mfem_error("Abstract method MultiBlockSolver::AddRHSFunction!\n"); }
   virtual void AddRHSFunction(const Vector &F)
   { mfem_error("Abstract method MultiBlockSolver::AddRHSFunction!\n"); }

   virtual void Assemble() = 0;
   virtual void AssembleRHS() = 0;
   virtual void AssembleOperator() = 0;
   // For bilinear case.
   // system-specific.
   virtual void AssembleInterfaceMatrices() = 0;

   // Global ROM operator Loading.
   virtual void SaveROMOperator(const std::string filename)
   { rom_handler->SaveOperator(filename); }
   virtual void LoadROMOperatorFromFile(const std::string filename)
   { rom_handler->LoadOperatorFromFile(filename); }

   // Component-wise assembly
   void GetComponentFESpaces(Array<FiniteElementSpace *> &comp_fes);
   // virtual void AllocateROMLinElems();

   void BuildROMLinElems();
   virtual void BuildCompROMLinElems() = 0;
   virtual void BuildBdrROMLinElems() = 0;
   // TODO(kevin): part of this can be transferred to InterfaceForm.
   virtual void BuildItfaceROMLinElems() = 0;

   void SaveROMLinElems(const std::string &filename)
   { assert(rom_elems); rom_elems->Save(filename); }

   void LoadROMLinElems(const std::string &filename)
   { assert(rom_elems); rom_elems->Load(filename); }

   void AssembleROMMat();

   virtual bool Solve(SampleGenerator *sample_generator = NULL) = 0;

   virtual void InitVisualization(const std::string& output_dir = "");
   virtual void InitUnifiedParaview(const std::string &file_prefix);
   virtual void InitIndividualParaview(const std::string &file_prefix);
   /* time-independent visualization */
   virtual void SaveVisualization();
   /* time-dependent visualization */
   virtual void SaveVisualization(const int step, const double time);

   void SaveSolution(std::string filename = "");
   void SaveSolutionWithTime(const std::string filename, const int step, const double time);
   void LoadSolution(const std::string &filename);
   void LoadSolutionWithTime(const std::string &filename, int &step, double &time);
   void CopySolution(BlockVector *input_sol);

   virtual void AllocateROMNlinElems()
   { mfem_error("Abstract method MultiBlockSolver::AllocateROMNlinElems!\n"); }
   virtual void BuildROMTensorElems()
   { mfem_error("Abstract method MultiBlockSolver::BuildROMTensorElems!\n"); }
   virtual void TrainROMEQPElems(SampleGenerator *sample_generator)
   { mfem_error("Abstract method MultiBlockSolver::TrainROMEQPElems!\n"); }
   virtual void SaveROMNlinElems(const std::string &input_prefix)
   { mfem_error("Abstract method MultiBlockSolver::SaveROMNlinElems!\n"); }
   virtual void LoadROMNlinElems(const std::string &input_prefix)
   { mfem_error("Abstract method MultiBlockSolver::LoadROMNlinElems!\n"); }
   virtual void AssembleROMNlinOper()
   { mfem_error("Abstract method MultiBlockSolver::AssembleROMNlinOper!\n"); }

   virtual void InitROMHandler();
   void GetBasisTags(std::vector<BasisTag> &basis_tags);

   virtual BlockVector* PrepareSnapshots(std::vector<BasisTag> &basis_tags);
   void SaveSnapshots(SampleGenerator *sample_generator);

   virtual void LoadReducedBasis() { rom_handler->LoadReducedBasis(); }
   virtual void ProjectOperatorOnReducedBasis() = 0;
   virtual void ProjectRHSOnReducedBasis();
   virtual void SolveROM();
   virtual void SaveBasisVisualization()
   { rom_handler->SaveBasisVisualization(fes, var_names); }

   virtual void SetParameterizedProblem(ParameterizedProblem *problem);

   void ComputeSubdomainErrorAndNorm(GridFunction *fom_sol, GridFunction *rom_sol, double &error, double &norm);
   void ComputeRelativeError(Array<GridFunction *> fom_sols, Array<GridFunction *> rom_sols, Vector &error);
   void CompareSolution(BlockVector &test_U, Vector &error);

   virtual void SaveEQPCoords(const std::string &filename) {}

protected:
   virtual void AssembleROMMat(BlockMatrix &romMat);
   virtual void AssembleROMMat(SparseBlockMatrixBuilder &builder);
};

#endif
