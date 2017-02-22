//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Bryce Adelstein-Lelbach
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2017 Zahra Khatami 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>

#include <hpx/hpx_init.hpp>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/util/high_resolution_timer.hpp>
#include <boost/range.hpp>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <typeinfo>
#include <iterator>
#include <hpx/parallel/executors/sequential_execution_wrapper.hpp>
#include <hpx/parallel/seq_or_par.hpp>
#include <hpx/parallel/chunk_size_determination.hpp>
#include <hpx/parallel/prefetching_distance_determination.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>
#include <hpx/parallel/executors/adaptive_chunk_size.hpp>

#include <boost/range.hpp>

///////////////////////////////////////////////////////////////////////////////
// Command-line variables
bool header = true; // print csv heading
double k = 0.5;     // heat transfer coefficient
double dt = 1.;     // time step
double dx = 1.;     // grid spacing

inline std::size_t idx(std::size_t i, int dir, std::size_t size)
{
    if(i == 0 && dir == -1)
        return size-1;
    if(i == size-1 && dir == +1)
        return 0;

    HPX_ASSERT((i + dir) < size);

    return i + dir;
}

///////////////////////////////////////////////////////////////////////////////
//[stepper_2
struct stepper
{
    // Our partition type
    typedef hpx::shared_future<double> partition;

    // Our data for one time step
    typedef std::vector<partition> space;

    // Our operator
    static double heat(double left, double middle, double right)
    {
        return middle + (k*dt/(dx*dx)) * (left - 2*middle + right);
    }

    // do all the work on 'nx' data points for 'nt' time steps
    hpx::future<space> do_work(std::size_t nx, std::size_t nt)
    {
        using hpx::dataflow;
        using hpx::util::unwrapped;

        // U[t][i] is the state of position i at time t.
        std::vector<space> U(2);
        for (space& s : U)
            s.resize(nx);

        // Initial conditions: f(0, i) = i
        for (std::size_t i = 0; i != nx; ++i)
            U[0][i] = hpx::make_ready_future(double(i));

        auto Op = unwrapped(&stepper::heat);

        // Actual time step loop
        for (std::size_t t = 0; t != nt; ++t)
        {
            space const& current = U[t % 2];
            space& next = U[(t + 1) % 2];

            // WHEN U[t][i-1], U[t][i], and U[t][i+1] have been computed, THEN we
            // can compute U[t+1][i]
            for (std::size_t i = 0; i != nx; ++i)
            {
                next[i] = dataflow(
                        hpx::launch::async, Op,
                        current[idx(i, -1, nx)], current[i], current[idx(i, +1, nx)]
                    );
            }
        }

        // Now the asynchronous computation is running; the above for-loop does not
        // wait on anything. There is no implicit waiting at the end of each timestep;
        // the computation of each U[t][i] will begin when as soon as its dependencies
        // are ready and hardware is available.

        // Return the solution at time-step 'nt'.
        return hpx::when_all(U[nt % 2]);
    }
};


///////////////////////////////////////////////////////////////////////////////

int hpx_main(boost::program_options::variables_map& vm)
{    
    // Applying runtime options
    std::uint64_t nx = vm["nx"].as<std::uint64_t>();   // Number of grid points.
    std::uint64_t nt = vm["nt"].as<std::uint64_t>();   // Number of steps.

    if (vm.count("no-header"))
        header = false;
    
    // Initialization and lambda function for implementing [2] and [3]
    std::vector<double> stnc(nx, 1.0);
    auto time_range = boost::irange(0, 45);

    auto f = [&stnc, nx](std::size_t t) {
        for(int i = 0; i < 100; i++) {  // FIXME (clang bug) should be equal to nx          
            if(i == 0) {
                stnc[i] = stnc[i] + (k*dt/(dx*dx)) * (- 2*stnc[i] + stnc[i + 1]);
            }
            else if(i == nx - 1) {
                stnc[i] = stnc[i] + (k*dt/(dx*dx)) * (stnc[i - 1] - 2*stnc[i]);
            }
            else {
                stnc[i] = stnc[i] + (k*dt/(dx*dx)) * (stnc[i - 1] - 2*stnc[i] + stnc[i + 1]);
            }
        }        
    };
    
    ////////////////////////////////////////////////////////////////////////
    // [1] Original code implemantion with HPX
    // Create the stepper object
    stepper step;
    // Measure execution time.
    std::size_t t_origin = hpx::util::high_resolution_clock::now();

    // Execute nt time steps on nx grid points.
    hpx::future<stepper::space> result = step.do_work(nx, nt);
    stepper::space solution = result.get();
    hpx::wait_all(solution);

    std::size_t elapsed_origin = hpx::util::high_resolution_clock::now() - t_origin;

    ////////////////////////////////////////////////////////////////////////
    // [2] Efficient chunk size
    std::size_t t_chunk = hpx::util::high_resolution_clock::now();

    hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::adaptive_chunk_size()), 
        time_range.begin(), time_range.end(), f);

    std::size_t elapsed_chunk = hpx::util::high_resolution_clock::now() - t_chunk;

    ////////////////////////////////////////////////////////////////////////
    // [3] Prefetching:
    std::size_t prefetch_distance_factor = 2;
    auto policy = hpx::parallel::par;
 
    std::size_t t_prefetch = hpx::util::high_resolution_clock::now();

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, prefetch_distance_factor, stnc), 
        time_range.begin(), time_range.end(), f);

    std::size_t elapsed_prefetch = hpx::util::high_resolution_clock::now() - t_prefetch;

    ////////////////////////////////////////////////////////////////////////
    // Printing results
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "time_origin = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_origin;
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "\ntime_chunk = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_chunk;
    std::cout << std::left << std::setw(6) << std::setfill(' ') << "\ntime_prefetching = ";
    std::cout << std::left << std::setw(6) << std::setfill(' ') << elapsed_prefetch << std::endl;
  
    return hpx::finalize();
}


int main(int argc, char* argv[])
{
    using namespace boost::program_options;

    options_description desc_commandline;
    desc_commandline.add_options()
        ("results", "print generated results (default: false)")
        ("nx", value<std::uint64_t>()->default_value(100),
         "Local x dimension")
        ("nt", value<std::uint64_t>()->default_value(45),
         "Number of time steps")
        ("k", value<double>(&k)->default_value(0.5),
         "Heat transfer coefficient (default: 0.5)")
        ("dt", value<double>(&dt)->default_value(1.0),
         "Timestep unit (default: 1.0[s])")
        ("dx", value<double>(&dx)->default_value(1.0),
         "Local x dimension")
        ( "no-header", "do not print out the csv header row")
    ;
    
    return hpx::init(desc_commandline, argc, argv);
}