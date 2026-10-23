// Copyright 2023 Lawrence Livermore National Security, LLC. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef ETC_HPP
#define ETC_HPP

#include "mfem.hpp"
#include "input_parser.hpp"

using namespace mfem;
using namespace std;

double UniformRandom();
int UniformRandom(const int &min, const int &max);

template <typename T>
inline void DeletePointers(Array<T*> &ptr_array)
{ for (int k = 0; k < ptr_array.Size(); k++) delete ptr_array[k]; }

template <typename T>
inline void DeletePointers(Array2D<T*> &ptr_array)
{
    for (int i = 0; i < ptr_array.NumRows(); i++)
        for (int j = 0; j < ptr_array.NumCols(); j++)
            delete ptr_array(i,j);
}

bool FileExists(const std::string& name);

/// Lightweight, non-collective timing output for parallel scaling studies.
bool ParallelProfileEnabled();
void PrintParallelProfile(const std::string &phase, double seconds,
                          long work_items = -1);

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

class TimeProfiler
{
private:
    std::string title = "";

    Array<StopWatch *> timers;

    Array<int> calls;
    Array<bool> starts;
    std::vector<std::string> names;
    std::unordered_map<std::string, int> indices;

    const MPI_Comm comm;
    int rank;

public:
    TimeProfiler(const std::string &title_, MPI_Comm comm_=MPI_COMM_WORLD)
        : title(title_), comm(comm_)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        calls.SetSize(0);
        names.clear();
        indices.clear();
    }

    virtual ~TimeProfiler()
    {
        std::string option = "time_profile/" + title;
        if (config.GetOption<bool>(option, false))
            Print(title);
        DeletePointers(timers);
    }

    void SetTitle(const std::string &title_)
    { title = title_; }

    void Start(const std::string &name)
    {
        if (!indices.count(name))
        {
            indices[name] = timers.Size();
            timers.Append(new StopWatch);
            names.push_back(name);
            calls.Append(0);
            starts.Append(false);
        }

        assert(indices.count(name));

        int idx = indices[name];
        if (starts[idx])
        {
            printf("TimeProfiler: %s is already started!\n", name.c_str());
            assert(!starts[idx]);
        }

        timers[idx]->Start();
        starts[idx] = true;
    }

    void Stop(const std::string &name)
    {
        assert(indices.count(name));
        
        int idx = indices[name];
        if (!starts[idx])
        {
            printf("TimeProfiler: %s is not started!\n", name.c_str());
            assert(starts[idx]);
        }

        timers[idx]->Stop();
        calls[idx] += 1;
        starts[idx] = false;
    }

    void Print(const std::string &title_, const bool compute_sum=false)
    {
        int nfunc = timers.Size();
        for (int k = 0; k < nfunc; k++)
            assert(!starts[k]);

        Array<double> times(nfunc);
        for (int k = 0; k < nfunc; k++)
            times[k] = timers[k]->RealTime();
        double total = times.Sum();

        MPI_Reduce(MPI_IN_PLACE, times.GetData(), nfunc, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(MPI_IN_PLACE, calls.GetData(), nfunc, MPI_INT, MPI_SUM, 0, comm);

        if (rank == 0)
        {
            printf("%s", (title_ + "\n").c_str());

            std::string line = std::string(100, '=');
            line += "\n";
            printf("%s", line.c_str());
            printf("%20s\t%20s\t%20s\t%20s\n", "Function", "Total time", "Calls", "Time per call");
            for (int k = 0; k < nfunc; k++)
            {
                printf("%20s\t%20.5e\t%20d\t%20.5e\n", names[k].c_str(), times[k], calls[k], times[k] / calls[k]);
            }
            if (compute_sum)
                printf("%20s\t%20.5e\t%20d\t%20.5e\n", "Total time", total, 1, total);
            printf("%s", line.c_str());
        }
    }
};

#endif
