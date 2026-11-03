// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include "main_workflow.hpp"
#include <cmath>

using namespace std;
using namespace mfem;

static const double threshold = 1.0e-14;
static const double stokes_threshold = 1.0e-11;
static const double linelast_threshold = 1.0e-13;

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound)
{
   SUCCEED();
}

TEST(Poisson_Workflow, SubmeshTest)
{
   config = InputParser("inputs/test.base.yml");

   config.dict_["visualization"]["enabled"] = true;

   config.dict_["model_reduction"]["rom_handler_type"] = "mfem";
   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold);

   return;
}

TEST(Poisson_Workflow, ComponentWiseTest)
{
   config = InputParser("inputs/test.component.yml");

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["mesh"]["type"] = "component-wise";
   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold);

   return;
}

TEST(Poisson_Workflow, ComponentWiseWithDirectSolve)
{
   config = InputParser("inputs/test.component.yml");
   config.dict_["solver"]["direct_solve"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "spd";
   config.dict_["model_reduction"]["hypre_matrix"] = "yes";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["mesh"]["type"] = "component-wise";
   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold);

   return;
}

TEST(Stokes_Workflow, SubmeshTest)
{
   config = InputParser("inputs/stokes.base.yml");

   config.dict_["visualization"]["enabled"] = true;

   config.dict_["model_reduction"]["rom_handler_type"] = "mfem";
   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   config.dict_["main"]["use_rom"] = false;
   GenerateSamples(MPI_COMM_WORLD);
   config.dict_["main"]["use_rom"] = true;

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(Stokes_Workflow, MFEMGlobalSeparateTest)
{
   config = InputParser("inputs/stokes.base.yml");

   config.dict_["visualization"]["enabled"] = true;

   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";
   config.dict_["model_reduction"]["linear_solver_type"] = "minres";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   config.dict_["main"]["use_rom"] = false;
   GenerateSamples(MPI_COMM_WORLD);
   config.dict_["main"]["use_rom"] = true;

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(Stokes_Workflow, ComponentWiseTest)
{
   config = InputParser("inputs/stokes.component.yml");

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   config.dict_["main"]["use_rom"] = false;
   GenerateSamples(MPI_COMM_WORLD);
   config.dict_["main"]["use_rom"] = true;

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(Stokes_Workflow, ComponentWiseWithDirectSolve)
{
   config = InputParser("inputs/stokes.component.yml");
   config.dict_["solver"]["direct_solve"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "sid";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(Stokes_Workflow, ComponentSeparateVariable)
{
   config = InputParser("inputs/stokes.component.yml");
   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["solver"]["direct_solve"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(Stokes_Workflow, ROM_OrderingByVariable)
{
   config = InputParser("inputs/stokes.component.yml");

   config.dict_["model_reduction"]["ordering"] = "variable";

   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["solver"]["direct_solve"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, SubmeshTest)
{
   config = InputParser("inputs/steady_ns.base.yml");

   config.dict_["mesh"]["uniform_refinement"] = 2;
   config.dict_["discretization"]["order"] = 2;
   config.dict_["visualization"]["enabled"] = true;

   config.dict_["model_reduction"]["rom_handler_type"] = "mfem";
   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, MFEMGlobalUniversalTest)
{
   config = InputParser("inputs/steady_ns.base.yml");

   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["mesh"]["uniform_refinement"] = 2;
   config.dict_["discretization"]["order"] = 2;
   config.dict_["visualization"]["enabled"] = true;

   config.dict_["model_reduction"]["rom_handler_type"] = "mfem";
   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, ComponentWiseTest)
{
   config = InputParser("inputs/steady_ns.component.yml");

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, ComponentWiseWithDirectSolve)
{
   config = InputParser("inputs/steady_ns.component.yml");
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, ComponentSeparateVariable)
{
   config = InputParser("inputs/steady_ns.component.yml");
   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD, "test_output.h5");

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, ComponentSeparateVariable_EQP)
{
   config = InputParser("inputs/steady_ns.component.yml");
   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";
   config.dict_["model_reduction"]["nonlinear_handling"] = "eqp";
   config.dict_["model_reduction"]["eqp"]["relative_tolerance"] = 1.0e-11;
   config.dict_["model_reduction"]["eqp"]["precompute"] = true;

   printf("\nSample Generation \n\n");
   
   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_eqp";
   TrainEQP(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD, "test_output.h5");

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(SteadyNS_Workflow, ROM_OrderingByVariable)
{
   config = InputParser("inputs/steady_ns.component.yml");

   config.dict_["model_reduction"]["ordering"] = "variable";

   config.dict_["model_reduction"]["separate_variable_basis"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "us";
   config.dict_["model_reduction"]["nonlinear_handling"] = "eqp";
   config.dict_["model_reduction"]["eqp"]["relative_tolerance"] = 1.0e-11;
   config.dict_["model_reduction"]["eqp"]["precompute"] = true;

   printf("\nSample Generation \n\n");
   
   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_eqp";
   TrainEQP(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD, "test_output.h5");

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < stokes_threshold);

   return;
}

TEST(LinElast_Workflow, SubmeshTest)
{
   config = InputParser("inputs/linelast.base.yml");

   config.dict_["model_reduction"]["rom_handler_type"] = "mfem";

   // Test save/loadSolution as well.
   config.dict_["save_solution"]["enabled"] = true;
   config.dict_["model_reduction"]["compare_solution"]["load_solution"] = true;
   config.dict_["model_reduction"]["compare_solution"]["fom_solution_file"] = "./sample1_solution.h5";

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);
   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["save_solution"]["enabled"] = false;
   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < linelast_threshold);

   return;
} 

TEST(LinElast_Workflow, ComponentWiseWithDirectSolve)
{
   config = InputParser("inputs/linelast.component.yml");
   config.dict_["solver"]["direct_solve"] = true;
   config.dict_["model_reduction"]["linear_solver_type"] = "direct";
   config.dict_["model_reduction"]["linear_system_type"] = "spd";

   printf("\nSample Generation \n\n");

    config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["mesh"]["type"] = "component-wise";
   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold); 

   return;
}

TEST(LinElast_Workflow, ComponentWiseTest)
{
   config = InputParser("inputs/linelast.component.yml");

   printf("\nSample Generation \n\n");

   config.dict_["main"]["mode"] = "sample_generation";
   GenerateSamples(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   printf("\nBuild ROM \n\n");

   config.dict_["mesh"]["type"] = "component-wise";
   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold);

   return;
}

TEST(AdvDiff_Workflow, SubmeshTest)
{
   config = InputParser("inputs/advdiff.base.yml");

   config.dict_["model_reduction"]["visualization"]["enabled"] = true;
   config.dict_["model_reduction"]["visualization"]["prefix"] = "basis_paraview";

   config.dict_["main"]["mode"] = "sample_generation";
   config.dict_["main"]["use_rom"] = false;
   GenerateSamples(MPI_COMM_WORLD);
   config.dict_["main"]["use_rom"] = true;

   config.dict_["main"]["mode"] = "train_rom";
   TrainROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "build_rom";
   BuildROM(MPI_COMM_WORLD);

   config.dict_["main"]["mode"] = "single_run";
   double error = SingleRun(MPI_COMM_WORLD);

   // This reproductive case must have a very small error at the level of finite-precision.
   printf("Error: %.15E\n", error);
   EXPECT_TRUE(error < threshold);

   return;
}

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   MPI_Init(&argc, &argv);
   int result = RUN_ALL_TESTS();
   MPI_Finalize();
   return result;
}
