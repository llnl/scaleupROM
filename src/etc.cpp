// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include "etc.hpp"
// #include <stdlib.h>
#include <random>

using namespace std;

static std::random_device rd;  // Will be used to obtain a seed for the random number engine
static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd(){}
static std::uniform_real_distribution<> dis(0.0, 1.0);

double UniformRandom()
{
   return dis(gen);
}

int UniformRandom(const int &min, const int &max)
{
   std::uniform_int_distribution<int> dis_int(min, max);
   return dis_int(gen);
}

bool FileExists(const std::string& name)
{
   std::ifstream f(name.c_str());
   return f.good();
   // ifstream f will be closed upon the end of the function.
}

bool ParallelProfileEnabled()
{
   return config.GetOption<bool>("time_profile/parallel_scaling", false);
}

void PrintParallelProfile(const std::string &phase, double seconds,
                          long work_items)
{
   if (!ParallelProfileEnabled()) return;
   int rank, nranks;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   printf("PARALLEL_PROFILE,%d,%d,%s,%.9g,%ld\n", rank, nranks,
          phase.c_str(), seconds, work_items);
   fflush(stdout);
}
