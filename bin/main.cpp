#include "main_workflow.hpp"
#include "etc.hpp"
#include <chrono>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   using Clock = std::chrono::steady_clock;
   const auto process_start = Clock::now();
   // 1. Initialize MPI.
    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   const auto mpi_initialized = Clock::now();

   const auto config_start = Clock::now();

   const char *input_file = "test.input";
   const char *output_file = "";
   const char *forced_input = "";
   OptionsParser args(argc, argv);
   args.AddOption(&input_file, "-i", "--input", "Input file to use.");
   args.AddOption(&output_file, "-o", "--output", "For single run, save comparison result to output file.");
   args.AddOption(&forced_input, "-f", "--forced-input",
      "Input options to overwrite. In the format of 'key1=value1:key2=value2:...'");
   args.ParseCheck();
   config = InputParser(input_file, forced_input);
   PrintParallelProfile(
      "app.mpi_init",
      std::chrono::duration<double>(mpi_initialized - process_start).count());
   PrintParallelProfile(
      "app.options_and_config",
      std::chrono::duration<double>(Clock::now() - config_start).count());

   std::string mode = config.GetOption<std::string>("main/mode", "run_example");
   const auto workflow_start = Clock::now();

   if (mode == "run_example")
   { if (rank == 0) RunExample(); }
   else 
   {
      if (mode == "sample_generation")    GenerateSamples(MPI_COMM_WORLD);
      else if (mode == "build_rom")       BuildROM(MPI_COMM_WORLD);
      else if (mode == "train_rom")       TrainROM(MPI_COMM_WORLD);
      else if (mode == "aux_train_rom")   AuxiliaryTrainROM(MPI_COMM_WORLD);
      else if (mode == "train_eqp")       TrainEQP(MPI_COMM_WORLD);
      else if (mode == "single_run")      SingleRun(MPI_COMM_WORLD, output_file);
      else if (mode == "print_eqp")       PrintEQPCoords(MPI_COMM_WORLD);
      else
      {
         if (rank == 0) printf("Unknown mode %s!\n", mode.c_str());
         exit(-1);
      }
   }

   PrintParallelProfile(
      "app.workflow",
      std::chrono::duration<double>(Clock::now() - workflow_start).count());
   const bool profile_enabled = ParallelProfileEnabled();
   const auto finalize_start = Clock::now();
   MPI_Finalize();
   const auto process_end = Clock::now();
   if (profile_enabled)
   {
      printf("PARALLEL_PROFILE,%d,%d,app.mpi_finalize,%.9g,-1\n",
             rank, num_procs,
             std::chrono::duration<double>(process_end - finalize_start).count());
      printf("PARALLEL_PROFILE,%d,%d,app.process_total,%.9g,-1\n",
             rank, num_procs,
             std::chrono::duration<double>(process_end - process_start).count());
      fflush(stdout);
   }
}
