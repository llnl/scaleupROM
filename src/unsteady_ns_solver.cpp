// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "unsteady_ns_solver.hpp"
#include "etc.hpp"

using namespace std;
using namespace mfem;

/*
   UnsteadyNSSolver
*/

UnsteadyNSSolver::UnsteadyNSSolver(TopologyHandler *input_topol_handler)
   : SteadyNSSolver(input_topol_handler), timer("UnsteadyNSSolver")
{
   nonlinear_mode = true;

   /* unsteady ns solver only uses LaxFridriechs flux solver */
   if (oper_type != OperType::LF)
      mfem_error("UnsteadyNSSolver only support Lax-Fridriechs flux integrator!\n");

   nt = config.GetRequiredOption<int>("time-integration/number_of_timesteps");
   dt = config.GetRequiredOption<double>("time-integration/timestep_size");
   time_order = config.GetOption<int>("time-integration/bdf_order", 1);
   report_interval = config.GetOption<int>("time-integration/report_interval", 0);
   metrics_file = config.GetOption<std::string>("metrics/filename", "metrics.h5");

   if (save_sol)
      restart_interval = config.GetOption<int>("save_solution/restart_interval", 0);

   /* optional min/max bounds checked at every timestep. an empty Vector disables the check. */
   Vector no_limit;
   u_limit = config.GetOption<Vector>("solution_limit/velocity", no_limit);
   p_limit = config.GetOption<Vector>("solution_limit/pressure", no_limit);

   if ((u_limit.Size() != 0) && (u_limit.Size() != 2))
      mfem_error("UnsteadyNSSolver: solution_limit/velocity must be a 2-element [min, max] list!\n");
   if ((p_limit.Size() != 0) && (p_limit.Size() != 2))
      mfem_error("UnsteadyNSSolver: solution_limit/pressure must be a 2-element [min, max] list!\n");
   if ((u_limit.Size() == 2) && (u_limit[0] > u_limit[1]))
      mfem_error("UnsteadyNSSolver: solution_limit/velocity must satisfy min <= max!\n");
   if ((p_limit.Size() == 2) && (p_limit[0] > p_limit[1]))
      mfem_error("UnsteadyNSSolver: solution_limit/pressure must satisfy min <= max!\n");

   if (time_order != 1)
      mfem_error("UnsteadyNSSolver supports only first-order time integration for now!\n");

   if (use_rom && !separate_variable_basis)
      mfem_error("UnsteadyNSSolver does not allow unified basis for all variables!\n");
}

UnsteadyNSSolver::~UnsteadyNSSolver()
{
   DeletePointers(mass);
   delete massMat;
   delete uu;
   delete U_step;
   delete RHS_step;
   delete U_stepview;
   delete RHS_stepview;
   delete Hop;
   delete rom_mass;

   delete u_ic;
   delete p_ic;
}

void UnsteadyNSSolver::BuildDomainOperators()
{
   SteadyNSSolver::BuildDomainOperators();

   mass.SetSize(numSub);
   for (int m = 0; m < numSub; m++)
   {
      mass[m] = new BilinearForm(ufes[m]);
      mass[m]->AddDomainIntegrator(new VectorMassIntegrator);
   }
}

void UnsteadyNSSolver::AssembleOperator()
{
   SteadyNSSolver::AssembleOperator();

   for (int m = 0; m < numSub; m++)
   {
      mass[m]->Assemble();
      mass[m]->Finalize();
   }

   massMat = new BlockMatrix(u_offsets);
   for (int i = 0; i < numSub; i++)
      massMat->SetBlock(i, i, &(mass[i]->SpMat()));
}

bool UnsteadyNSSolver::Solve(SampleGenerator *sample_generator)
{
   bool converged = true;
   int initial_step = 0;
   double time = 0.0;

   std::string restart_file, file_fmt;
   file_fmt = "%s/%s_%08d.h5";

   timer.Start("Solve/init");

   SetupInitialCondition(initial_step, time);

   /* physical time at which this run starts, nonzero when restarting */
   const double initial_time = time;

   int sample_interval = config.GetOption<int>("sample_generation/time-integration/sample_interval", 0);
   int bootstrap = config.GetOption<int>("sample_generation/time-integration/bootstrap", 0);

   InitializeTimeIntegration();

   SortByVariables(*U, *U_step);

   SaveVisualization(0, time);

   timer.Stop("Solve/init");

   double cfl = 0.0;
   for (int step = initial_step; step < nt; step++)
   {
      Step(time, step);

      timer.Start("Solve/save_sol");

      cfl = ComputeCFL(dt);
      SanityCheck(step, time - initial_time);
      if (report_interval &&
          ((step+1) % report_interval) == 0)
         printf("Time step: %05d, CFL: %.3e\n", step+1, cfl);

      if (visual.time_interval &&
          ((step+1) % visual.time_interval) == 0)
         SaveVisualization(step+1, time);

      if (save_sol && restart_interval &&
          (((step+1) % restart_interval) == 0))
      {
         restart_file = string_format(file_fmt, sol_dir.c_str(), sol_prefix.c_str(), step+1);
         SaveSolutionWithTime(restart_file, step+1, time);
      }

      /* save solution if sample generator is provided */
      if (sample_generator && sample_interval &&
          ((step+1) > bootstrap) && (((step+1) % sample_interval) == 0))
         SaveSnapshots(sample_generator);

      timer.Stop("Solve/save_sol");
   }

   timer.Start("Solve/final");

   SortBySubdomains(*U_step, *U);

   timer.Stop("Solve/final");

   timer.Print("UnsteadyNSSolver::Solve", true);

   SaveMetrics(converged, time - initial_time);

   return converged;
}

void UnsteadyNSSolver::SanityCheck(const int step, const double simulation_time)
{
   if (isnan(U_step->Min()) || isnan(U_step->Max()))
   {
      printf("Step : %d\n", step);
      SaveMetrics(false, simulation_time);
      mfem_error("UnsteadyNSSolver: Solution blew up!!\n");
   }

   if ((u_limit.Size() == 0) && (p_limit.Size() == 0)) return;

   /* U_step is ordered by variable: block m is subdomain m's velocity,
      block numSub + m its pressure. */
   for (int m = 0; m < numSub; m++)
   {
      if (u_limit.Size() == 2)
      {
         Vector &u_m = U_step->GetBlock(m);

         /* |u| <= sqrt(dim) * max_d |u_d|, and a magnitude never falls below a non-positive
            bound, so these are the only cases where a node can be out of range. */
         if ((sqrt((double) vdim[0]) * u_m.Normlinf() > u_limit[1]) || (u_limit[0] > 0.0))
         {
            const int ndofs = ufes[m]->GetNDofs();
            /* ufes is Ordering::byNODES and DenseMatrix is column-major, so column d of this
               view is velocity component d and row n is node n's velocity vector. */
            DenseMatrix u_nodal(u_m.GetData(), ndofs, vdim[0]);
            Vector u_n;
            for (int n = 0; n < ndofs; n++)
            {
               u_nodal.GetRow(n, u_n);
               const double umag = u_n.Norml2();
               if ((umag >= u_limit[0]) && (umag <= u_limit[1])) continue;

               ReportSolutionCrash(step, simulation_time, m, 0, n, umag, u_limit);
            }
         }
      }

      if (p_limit.Size() == 2)
      {
         Vector &p_m = U_step->GetBlock(numSub + m);

         if ((p_m.Min() < p_limit[0]) || (p_m.Max() > p_limit[1]))
            for (int n = 0; n < p_m.Size(); n++)
            {
               if ((p_m(n) >= p_limit[0]) && (p_m(n) <= p_limit[1])) continue;

               ReportSolutionCrash(step, simulation_time, m, 1, n, p_m(n), p_limit);
            }
      }
   }
}

void UnsteadyNSSolver::ReportSolutionCrash(const int step, const double simulation_time, const int m,
                                           const int var, const int dof, const double val,
                                           const Vector &limit)
{
   /* Mesh::GetNodes projects the coordinates onto the given GridFunction's space. A dim-component
      space on fec[var] shares the scalar dof numbering of ufes[m]/pfes[m], so DofToVDof(dof, d)
      is coordinate d of the offending node. */
   FiniteElementSpace coord_fes(meshes[m], fec[var], dim);
   GridFunction coords(&coord_fes);
   meshes[m]->GetNodes(coords);

   /* every rank runs the identical time integration, so only rank 0 reports. */
   if (rank == 0)
   {
      printf("UnsteadyNSSolver: solution crashed at step %d (t = %.4e)!\n", step + 1, simulation_time);
      printf("  subdomain %d, %s = %.4e outside [%.4e, %.4e] at (", m,
             (var == 0) ? "|vel|" : "pres", val, limit[0], limit[1]);
      for (int d = 0; d < dim; d++)
         printf(" %.4e", coords(coord_fes.DofToVDof(dof, d)));
      printf(" )\n");
   }

   SaveMetrics(false, simulation_time);
   mfem_error("UnsteadyNSSolver: Solution went beyond solution_limit!\n");
}

void UnsteadyNSSolver::SaveMetrics(const bool converged, const double simulation_time)
{
   /* every rank runs the identical time integration, so only rank 0 writes the file. */
   if (rank != 0) return;

   printf("Saving the metrics file %s ...", metrics_file.c_str());

   hid_t file_id;
   herr_t errf = 0;
   file_id = H5Fcreate(metrics_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   assert(file_id >= 0);

   hdf5_utils::WriteAttribute(file_id, "converged", converged);
   hdf5_utils::WriteAttribute(file_id, "simulation_time", simulation_time);

   errf = H5Fclose(file_id);
   assert(errf >= 0);
   printf("Done!\n");
}

void UnsteadyNSSolver::InitializeTimeIntegration()
{
   assert(dt > 0.0);

   uu = new SparseMatrix(vblock_offsets[1]);
   SparseMatrix *tmp = massMat->CreateMonolithic();
   (*uu) += *tmp;
   (*uu) *= bd0 / dt;
   (*uu) += (*M); // add viscous flux operator
   uu->Finalize();
   delete tmp;

   // This should be run after AssembleOperator
   assert(systemOp);
   systemOp->SetBlock(0, 0, uu);
   systemOp->SetBlock(0, 1, Bt);
   systemOp->SetBlock(1, 0, B);

   StokesSolver::SetupMUMPSSolver(true);

   Hop = new BlockOperator(u_offsets);
   for (int m = 0; m < numSub; m++)
      Hop->SetDiagonalBlock(m, hs[m]);

   /* the routines above can be merged to AssembleOperator */

   offsets_byvar.SetSize(num_var * numSub + 1);
   offsets_byvar = 0;
   for (int k = 0; k < numSub; k++)
   {
      offsets_byvar[k+1] = u_offsets[k+1];
      offsets_byvar[k+1 + numSub] = p_offsets[k+1] + u_offsets.Last();
   }

   U_step = new BlockVector(offsets_byvar);
   RHS_step = new BlockVector(offsets_byvar);
   U_stepview = new BlockVector(U_step->GetData(), vblock_offsets);
   RHS_stepview = new BlockVector(RHS_step->GetData(), vblock_offsets);

   u1.SetSize(U_stepview->BlockSize(0));
   Cu1.SetSize(U_stepview->BlockSize(0));
}

void UnsteadyNSSolver::SetupInitialCondition(int &initial_step, double &time)
{
   bool use_restart = config.GetOption<bool>("solver/use_restart", false);
   std::string restart_file, file_fmt;
   file_fmt = "%s/%s_%08d.h5";

   if (use_restart)
   {
      int restart_timestep = config.GetOption<int>("solver/restart_timestep", -1);
      if (restart_timestep > 0)
      {
         std::string backup_file = string_format(file_fmt, sol_dir.c_str(), sol_prefix.c_str(), restart_timestep);
         restart_file = config.GetOption<std::string>("solver/restart_file", backup_file);
      }
      else
         restart_file = config.GetRequiredOption<std::string>("solver/restart_file");
      LoadSolutionWithTime(restart_file, initial_step, time);
   }
   else
   {
      if (u_ic && p_ic)
         for (int m = 0; m < numSub; m++)
         {
            vels[m]->ProjectCoefficient(*u_ic);
            ps[m]->ProjectCoefficient(*p_ic);
         }

      initial_step = 0;
      time = 0.0;
   }

   if ((!use_restart) && save_sol)
   {
      restart_file = string_format(file_fmt, sol_dir.c_str(), sol_prefix.c_str(), initial_step);
      SaveSolutionWithTime(restart_file, initial_step, time);
   }
}

void UnsteadyNSSolver::Step(double &time, int step)
{
   timer.Start("Step/copy_u");

   /* set time for forcing/boundary. At this point, time remains at the previous timestep. */
   SetTime(time);

   /* copy velocity */
   u1 = U_stepview->GetBlock(0);

   timer.Stop("Step/copy_u");
   timer.Start("Step/advect_dom");

   /* evaluate nonlinear advection at previous time step */
   Hop->Mult(u1, Cu1);

   timer.Stop("Step/advect_dom");
   timer.Start("Step/advect_itf");

   nl_itf->InterfaceAddMult(u1, Cu1);

   timer.Stop("Step/advect_itf");
   timer.Start("Step/rhs_sort");

   /* Base right-hand side for boundary conditions and forcing */
   SortByVariables(*RHS, *RHS_step);

   timer.Stop("Step/rhs_sort");
   timer.Start("Step/add_advect");

   /* Add nonlinear convection */
   RHS_stepview->GetBlock(0).Add(-ab1, Cu1);

   timer.Stop("Step/add_advect");
   timer.Start("Step/add_udot");

   /* Add time derivative term */
   // TODO: extend for high order bdf schemes
   massMat->AddMult(u1, RHS_stepview->GetBlock(0), -bd1 / dt);

   timer.Stop("Step/add_udot");
   timer.Start("Step/solve_system");

   /* Solve for the next step */
   mumps->Mult(*RHS_step, *U_step);

   timer.Stop("Step/solve_system");
   timer.Start("Step/pres_scale");

   /* remove pressure scalar if all dirichlet bc */
   if (!pres_dbc)
   {
      double p_const = U_stepview->GetBlock(1).Sum() / U_stepview->GetBlock(1).Size();

      U_stepview->GetBlock(1) -= p_const;
   }

   timer.Stop("Step/pres_scale");

   time += dt;
}

void UnsteadyNSSolver::SaveVisualization(const int step, const double time)
{
   /* copy to original solution variable */
   SortBySubdomains(*U_step, *U);

   MultiBlockSolver::SaveVisualization(step, time);
}

bool UnsteadyNSSolver::SetParameterizedProblem(ParameterizedProblem *problem)
{
   SteadyNSSolver::SetParameterizedProblem(problem);

   /* initialize solution with zero first */
   for (int m = 0; m < numSub; m++)
   {
      *vels[m] = 0.0;
      *ps[m] = 0.0;
   }

   Vector zero_vel(dim), zero_pres(1);
   zero_vel = 0.0; zero_pres = 0.0;

   if (u_ic) delete u_ic;
   if (p_ic) delete p_ic;

   /* set up initial condition */
   bool ic_setup = (problem->ic_ptr.Size() >= 2);
   if (!ic_setup)
   {
      u_ic = new VectorConstantCoefficient(zero_vel);
      p_ic = new VectorConstantCoefficient(zero_pres);
      // UnsteadyNSSolver does not fail in SetParameterizedProblem.
      return true;
   }

   if (problem->ic_ptr[0])
      u_ic = new VectorFunctionCoefficient(vdim[0], problem->ic_ptr[0]);
   else
      u_ic = new VectorConstantCoefficient(zero_vel);

   if (problem->ic_ptr[1])
      p_ic = new VectorFunctionCoefficient(1, problem->ic_ptr[1]);
   else
      p_ic = new VectorConstantCoefficient(zero_pres);
   // UnsteadyNSSolver does not fail in SetParameterizedProblem.
   return true;
}

BlockVector* UnsteadyNSSolver::PrepareSnapshots(std::vector<BasisTag> &basis_tags)
{
   /* copy to original solution variable */
   SortBySubdomains(*U_step, *U);

   return MultiBlockSolver::PrepareSnapshots(basis_tags);
}

double UnsteadyNSSolver::ComputeCFL(const double dt_)
{
   Vector ux, uy, uz;
   double cflx = 0.0;
   double cfly = 0.0;
   double cflz = 0.0;
   double cflm = 0.0;
   double cflmax = 0.0;

   for (int m = 0; m < numSub; m++)
   {
      GridFunction vel;
      vel.MakeRef(ufes[m], U_step->GetBlock(m), 0);

      for (int e = 0; e < ufes[m]->GetNE(); ++e)
      {
         const FiniteElement *fe = ufes[m]->GetFE(e);
         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                fe->GetOrder());
         ElementTransformation *tr = ufes[m]->GetElementTransformation(e);

         vel.GetValues(e, ir, ux, 1);
         vel.GetValues(e, ir, uy, 2);
         if (vdim[0] == 3)
            vel.GetValues(e, ir, uz, 3);

         double hmin = meshes[m]->GetElementSize(e, 1) /
                     (double) ufes[m]->GetElementOrder(0);

         for (int i = 0; i < ir.GetNPoints(); ++i)
         {
            const IntegrationPoint &ip = ir.IntPoint(i);
            tr->SetIntPoint(&ip);
            const DenseMatrix &invJ = tr->InverseJacobian();
            const double detJinv = 1.0 / tr->Jacobian().Det();

            cflx = fabs(dt_ * ux(i) / hmin);
            cfly = fabs(dt_ * uy(i) / hmin);
            if (vdim[0] == 3)
               cflz = fabs(dt_ * uz(i) / hmin);

            cflm = cflx + cfly + cflz;
            cflmax = fmax(cflmax, cflm);
         }  // for (int i = 0; i < ir.GetNPoints(); ++i)
      }  // for (int e = 0; e < ufes->GetNE(); ++e)
   }  // for (int m = 0; m < numSub; m++)

   double cflmax_global = cflmax;

   return cflmax_global;
}

void UnsteadyNSSolver::SetTime(const double time)
{
   /* set time for forcing coefficients */
   for (int k = 0; k < f_coeffs.Size(); k++)
      if (f_coeffs[k]) f_coeffs[k]->SetTime(time);

   /* set time for boundary conditions */
   for (int k = 0; k < ud_coeffs.Size(); k++)
      if (ud_coeffs[k]) ud_coeffs[k]->SetTime(time);
   for (int k = 0; k < sn_coeffs.Size(); k++)
      if (sn_coeffs[k]) sn_coeffs[k]->SetTime(time);
}

void UnsteadyNSSolver::InitROMHandler()
{
   SteadyNSSolver::InitROMHandler();

   if (rom_handler->GetOrdering() != ROMOrderBy::VARIABLE)
      mfem_error("UnsteadyNSSolver::InitROMHandler- unsteady NS solver only allows ROM ordering by variable!\n");
}

void UnsteadyNSSolver::BuildCompROMLinElems()
{
   SteadyNSSolver::BuildCompROMLinElems();

   if (!separate_variable_basis)
      mfem_error("UnsteadyNSSolver::BuildCompROMLinElems- unified basis is not supported!\n");

   /* We only need velocity mass matrix */
   const int num_comp = topol_handler->GetNumComponents();
   for (int c = 0; c < num_comp; c++)
   {
      const int fidx = c * num_var;
      Mesh *comp = topol_handler->GetComponentMesh(c);

      BilinearForm m_comp(comp_fes[fidx]);

      m_comp.AddDomainIntegrator(new VectorMassIntegrator);
      m_comp.Assemble();
      m_comp.Finalize();      

      SparseMatrix *m_mat = &(m_comp.SpMat());

      assert(rom_elems->mass[c]->nrows > 0);
      assert(rom_elems->mass[c]->ncols > 0);

      if (separate_variable_basis)
         (*rom_elems->mass[c])(0, 0) = rom_handler->ProjectToRefBasis(fidx, fidx, m_mat);
   }
}

void UnsteadyNSSolver::AssembleROMMat(BlockMatrix &romMat)
{
   SteadyNSSolver::AssembleROMMat(romMat);

   assert(rom_handler->GetOrdering() == ROMOrderBy::VARIABLE);

   rom_handler->GetBlockOffsets()->GetSubArray(0, numSub+1, rom_u_offsets);
   rom_mass = new BlockMatrix(rom_u_offsets);
   rom_mass->owns_blocks = true;

   for (int m = 0; m < numSub; m++)
   {
      int c_type = topol_handler->GetMeshType(m);

      Array<int> midx(1);
      midx[0] = rom_handler->GetBlockIndex(m, 0);

      MatrixBlocks mass_mat(1, 1);
      mass_mat(0, 0) = new SparseMatrix(*(*rom_elems->mass[c_type])(0,0));

      /* mass matrix itself for RHS */
      AddToBlockMatrix(midx, midx, mass_mat, *rom_mass);

      /* add mass matrices on LHS */
      *mass_mat(0, 0) *= bd0 / dt;
      AddToBlockMatrix(midx, midx, mass_mat, romMat);
   }  // for (int m = 0; m < numSub; m++)
}

void UnsteadyNSSolver::SolveROM()
{
   assert(rom_handler->GetOrdering() == ROMOrderBy::VARIABLE);

   timer.Start("SolveROM/init");

   int initial_step = 0;
   double time = 0.0;

   std::string restart_file, file_fmt;
   file_fmt = "%s/%s_%08d.h5";
   
   SetupInitialCondition(initial_step, time);   

   const Array<int> *rom_block_offsets = rom_handler->GetBlockOffsets();

   Array<int> rom_u_offsets(numSub + 1);
   Array<int> rom_p_offsets(numSub + 1);
   rom_block_offsets->GetSubArray(0, numSub + 1, rom_u_offsets);
   rom_block_offsets->GetSubArray(numSub, numSub + 1, rom_p_offsets);
   for (int k = numSub; k >= 0; k--)
      rom_p_offsets[k] -= rom_p_offsets[0];

   Array<int> rom_var_offsets(num_var + 1);
   rom_var_offsets[0] = 0;
   rom_var_offsets[1] = (*rom_block_offsets)[numSub];
   rom_var_offsets[2] = (*rom_block_offsets)[2 * numSub];

   BlockVector *reduced_sol = NULL;
   rom_handler->ProjectGlobalToDomainBasis(U, reduced_sol);
   BlockVector reduced_rhs(*rom_handler->GetReducedRHS());

   BlockVector *rsol_view = new BlockVector(reduced_sol->GetData(), rom_var_offsets);
   BlockVector *rrhs_view = new BlockVector(reduced_rhs.GetData(), rom_var_offsets);

   u1.SetSize(rsol_view->BlockSize(0));
   Cu1.SetSize(rsol_view->BlockSize(0));

   Hop = new BlockOperator(rom_u_offsets);
   for (int m = 0; m < numSub; m++)
      Hop->SetDiagonalBlock(m, subdomain_eqps[m]);

   int pN = -1;
   BlockVector rom_ones(rom_p_offsets);
   rom_ones = 0.0;
   if (!pres_dbc)
   {
      pN = p_offsets.Last();
      BlockVector fom_ones(p_offsets);
      fom_ones = 1.0;
      for (int m = 0; m < numSub; m++)
      {
         int idx = rom_handler->GetBlockIndex(m, 1);
         rom_handler->ProjectToDomainBasis(idx, fom_ones.GetBlock(m), rom_ones.GetBlock(m));
      }
   }

   timer.Stop("SolveROM/init");

   for (int step = initial_step; step < nt; step++)
   {
      timer.Start("SolveROM/set_time");

      /* set time for forcing/boundary. At this point, time remains at the previous timestep. */
      SetTime(time);

      timer.Stop("SolveROM/set_time");
      timer.Start("SolveROM/copy_u");

      /* copy velocity */
      u1 = rsol_view->GetBlock(0);

      timer.Stop("SolveROM/copy_u");
      timer.Start("SolveROM/advect_dom");

      /* evaluate nonlinear advection at previous time step */
      Hop->Mult(u1, Cu1);

      timer.Stop("SolveROM/advect_dom");
      timer.Start("SolveROM/advect_itf");

      itf_eqp->InterfaceAddMult(u1, Cu1);

      timer.Stop("SolveROM/advect_itf");
      timer.Start("SolveROM/get_rhs");

      /* Base right-hand side for boundary conditions and forcing */
      reduced_rhs = *(rom_handler->GetReducedRHS());

      timer.Stop("SolveROM/get_rhs");
      timer.Start("SolveROM/add_advect");

      /* Add nonlinear convection */
      rrhs_view->GetBlock(0).Add(-ab1, Cu1);

      timer.Stop("SolveROM/add_advect");
      timer.Start("SolveROM/add_udot");

      /* Add time derivative term */
      // TODO: extend for high order bdf schemes
      rom_mass->AddMult(u1, rrhs_view->GetBlock(0), -bd1 / dt);

      timer.Stop("SolveROM/add_udot");
      timer.Start("SolveROM/solve_system");

      /* Solve for the next step */
      rom_handler->Solve(reduced_rhs, *reduced_sol);

      timer.Stop("SolveROM/solve_system");
      timer.Start("SolveROM/pres_scale");

      /* remove pressure scalar if all dirichlet bc */
      if (!pres_dbc)
      {
         double p_const = (rom_ones * (rsol_view->GetBlock(1))) / pN;

         rsol_view->GetBlock(1).Add(-p_const, rom_ones);
      }

      timer.Stop("SolveROM/pres_scale");

      time += dt;
   }

   timer.Start("SolveROM/liftup");

   rom_handler->LiftUpGlobal(*reduced_sol, *U);

   timer.Stop("SolveROM/liftup");

   timer.Print("UnsteadyNSSolver::SolveROM", true);

   delete rsol_view;
   delete reduced_sol;
   return;
}