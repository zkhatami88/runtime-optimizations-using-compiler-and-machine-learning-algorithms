//  Copyright (c) 2014 Hartmut Kaiser
//  Copyright (c) 2014 Bryce Adelstein-Lelbach
//  Copyright (c) 2014 Patricia Grubel
//  Copyright (c) 2017 Zahra Khatami 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
/*
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

    
//DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::chunk_size_determination({hpx::get_os_thread_count(), 3502, 2500, 301, std::size_t(std::distance(time_range.begin(), time_range.end())), 1})),  time_range.begin(), time_range.end(), f);

    std::size_t elapsed_chunk = hpx::util::high_resolution_clock::now() - t_chunk;

    ////////////////////////////////////////////////////////////////////////
    // [3] Prefetching:
    std::size_t prefetch_distance_factor = 2;
    auto policy = hpx::parallel::par;
 
    std::size_t t_prefetch = hpx::util::high_resolution_clock::now();

    
//DETERMING PREFETCHER DISTANCE BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, hpx::parallel::prefetching_distance_determination({hpx::get_os_thread_count(), 3502, 2500, 301, std::size_t(std::distance(time_range.begin(), time_range.end())), 1}), stnc), 
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
*/

// Copyright (c) 2016 Thomas Heller
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include "stencil.hpp"
#include "output.hpp"
#include "communicator.hpp"

#include <hpx/hpx_init.hpp>
#include <hpx/include/compute.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/parallel_algorithm.hpp>

#include <hpx/util/high_resolution_timer.hpp>

#include <hpx/parallel/executors/sequential_execution_wrapper.hpp>
#include <hpx/parallel/seq_or_par.hpp>
#include <hpx/parallel/chunk_size_determination.hpp>
#include <hpx/parallel/prefetching_distance_determination.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>
#include <hpx/parallel/executors/adaptive_chunk_size.hpp>

#include <array>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

typedef hpx::compute::host::block_allocator<double> allocator_type;
typedef hpx::compute::host::block_executor<> executor_type;
typedef hpx::compute::vector<double, allocator_type> data_type;
typedef row_iterator<data_type::iterator> iterator;

typedef std::vector<double> communication_type;

HPX_REGISTER_CHANNEL_DECLARATION(communication_type);
HPX_REGISTER_CHANNEL(communication_type, stencil_communication);

void worker_original(
    std::size_t rank, std::size_t num, std::size_t Nx, std::size_t Ny, std::size_t steps,
    std::string const& output_name)
{
    std::array<data_type, 2> U;

    auto numa_domains = hpx::compute::host::numa_domains();
    allocator_type alloc(numa_domains);

    U[0] = data_type(Nx * Ny, 0.0, alloc);
    U[1] = data_type(Nx * Ny, 0.0, alloc);

    init(U, Nx, Ny, rank, num);

    // Setup our communicator
    typedef communicator<std::vector<double>> communicator_type;
    communicator_type comm(rank, num);

    if (comm.has_neighbor(communicator_type::up))
    {
        // send initial value to our upper neighbor
        comm.set(communicator_type::up,
            std::vector<double>(U[0].begin(), U[0].begin() + Nx), 0);
    }
    if (comm.has_neighbor(communicator_type::down))
    {
            // send initial value to our neighbor below
        comm.set(communicator_type::down,
            std::vector<double>(U[0].end() - Nx, U[0].end()), 0);
    }

    if (rank == 0)
    {
        std::cout << "Running original-example using " << num << " Partitions\n";
    }

    executor_type executor(numa_domains);
    hpx::util::high_resolution_timer t;

    // Construct our column iterators. We want to begin with the second
    // row to avoid out of bound accesses.
    iterator curr(Nx, U[0].begin());
    iterator next(Nx, U[1].begin());

    auto policy = hpx::parallel::par.on(executor);
    for (std::size_t t = 0; t < steps; ++t)
    {
        // Update our upper boundary if we have an interior partition and an
        // upper neighbor
        if (comm.has_neighbor(communicator_type::up))
        {
            // Get the first row.
            auto result = next.middle;
            // retrieve the row which is 'up' from our first row.
            std::vector<double> up = comm.get(communicator_type::up, t).get();
            // Create a row iterator with that top boundary
            auto it = curr.top_boundary(up);

            // After getting our missing row, we can update our first row
            line_update(it, it + Nx, result);

            // Finally, we can send the updated first row for our neighbor
            // to consume in the next timestep. Don't send if we are on
            // the last timestep
            comm.set(communicator_type::up,
                std::vector<double>(result, result + Nx), t + 1);
        }

        // Update our interior spatial domain
        hpx::parallel::for_loop(policy,
            curr + 1, curr + Ny-1,
            // We need to advance the result by one row each iteration
            hpx::parallel::induction(next.middle + Nx, Nx),
            [Nx](iterator it, data_type::iterator result)
            {
                line_update(*it, *it + Nx, result);
            }
        );

        // Update our lower boundary if we have an interior partition and a
        // neighbor below
        if (comm.has_neighbor(communicator_type::down))
        {
            // Get the last row.
            auto result = next.middle + (Ny - 2) * Nx;
            // retrieve the row which is 'down' from our last row.
            std::vector<double> down = comm.get(communicator_type::down, t).get();
            // Create a row iterator with that bottom boundary
            auto it = (curr + Ny - 2).bottom_boundary(down);
            // After getting our missing row, we can update our last row
            line_update(it, it + Nx, result);

            // Finally, we can send the updated last row for our neighbor
            // to consume in the next timestep. Don't send if we are on
            // the last timestep
            comm.set(communicator_type::down,
                std::vector<double>(result, result + Nx), t + 1);
        }

        if (rank == 0)
            std::cout << "." << std::flush;

        std::swap(curr, next);
    }
    double elapsed = t.elapsed();

    if (rank == 0)
        std::cout << "\n";

    if (rank == 0)
    {
        double mlups = (((Nx - 2.) * (Ny * num - 2.) * steps) / 1e6)/ elapsed;
        std::cout << "MLUPS: " << mlups << "\n";
    }

    if (!output_name.empty())
        output(output_name + std::to_string(rank), U[0], Nx, Ny);
}

///////////////////////////////////////////////////////////////////////////////
// Using efficient chunk size

void worker_chunk(
    std::size_t rank, std::size_t num, std::size_t Nx, std::size_t Ny, std::size_t steps,
    std::string const& output_name)
{
    std::array<data_type, 2> U;

    auto numa_domains = hpx::compute::host::numa_domains();
    allocator_type alloc(numa_domains);

    U[0] = data_type(Nx * Ny, 0.0, alloc);
    U[1] = data_type(Nx * Ny, 0.0, alloc);

    init(U, Nx, Ny, rank, num);

    // Setup our communicator
    typedef communicator<std::vector<double>> communicator_type;
    communicator_type comm(rank, num);

    if (comm.has_neighbor(communicator_type::up))
    {
        // send initial value to our upper neighbor
        comm.set(communicator_type::up,
            std::vector<double>(U[0].begin(), U[0].begin() + Nx), 0);
    }
    if (comm.has_neighbor(communicator_type::down))
    {
            // send initial value to our neighbor below
        comm.set(communicator_type::down,
            std::vector<double>(U[0].end() - Nx, U[0].end()), 0);
    }

    if (rank == 0)
    {
        std::cout << "Running efficient-chunk-example using " << num << " Partitions\n";
    }

    executor_type executor(numa_domains);
    

    // Construct our column iterators. We want to begin with the second
    // row to avoid out of bound accesses.
    iterator curr(Nx, U[0].begin());
    iterator next(Nx, U[1].begin());

    auto policy = hpx::parallel::par.on(executor);
    auto f = [Nx](iterator it, data_type::iterator result) {
                line_update(*it, *it + Nx, result);
            };

    hpx::util::high_resolution_timer t;      
    for (int t = 0; t < steps; ++t)
    {
        // Update our upper boundary if we have an interior partition and an
        // upper neighbor
        if (comm.has_neighbor(communicator_type::up))
        {
            // Get the first row.
            auto result = next.middle;
            // retrieve the row which is 'up' from our first row.
            std::vector<double> up = comm.get(communicator_type::up, t).get();
            // Create a row iterator with that top boundary
            auto it = curr.top_boundary(up);

            // After getting our missing row, we can update our first row
            line_update(it, it + Nx, result);

            // Finally, we can send the updated first row for our neighbor
            // to consume in the next timestep. Don't send if we are on
            // the last timestep
            comm.set(communicator_type::up,
                std::vector<double>(result, result + Nx), t + 1);
        }

        auto begin = curr + 1;
        auto end = curr + Ny - 1;
        
        // Update our interior spatial domain
        
//DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_loop(policy.with(hpx::parallel::chunk_size_determination({hpx::get_os_thread_count(), 0, 0, 0, std::size_t(std::distance(begin, end)), 0})),  begin, end, hpx::parallel::induction(next.middle + Nx, Nx), f);

        // Update our lower boundary if we have an interior partition and a
        // neighbor below
        if (comm.has_neighbor(communicator_type::down))
        {
            // Get the last row.
            auto result = next.middle + (Ny - 2) * Nx;
            // retrieve the row which is 'down' from our last row.
            std::vector<double> down = comm.get(communicator_type::down, t).get();
            // Create a row iterator with that bottom boundary
            auto it = (curr + Ny - 2).bottom_boundary(down);
            // After getting our missing row, we can update our last row
            line_update(it, it + Nx, result);

            // Finally, we can send the updated last row for our neighbor
            // to consume in the next timestep. Don't send if we are on
            // the last timestep
            comm.set(communicator_type::down,
                std::vector<double>(result, result + Nx), t + 1);
        }

        if (rank == 0)
            std::cout << "." << std::flush;

        std::swap(curr, next);
    }
    double elapsed = t.elapsed();

    if (rank == 0)
        std::cout << "\n";

    if (rank == 0)
    {
        double mlups = (((Nx - 2.) * (Ny * num - 2.) * steps) / 1e6)/ elapsed;
        std::cout << "MLUPS_chunk: " << mlups << "\n";
    }

    if (!output_name.empty())
        output(output_name + std::to_string(rank), U[0], Nx, Ny);
}

///////////////////////////////////////////////////////////////////////////////
// Using efficient prefetcher distance
void line_update_prefetch(int begin, int end, std::vector<double>& U0, std::vector<double>& U1) {
    for(int i = begin; i <= end; i++) {
        U1[i] = 0.5 * (U0[i - 1] + U0[i + 1]);
    }
}

void worker_prefecth(
    std::size_t rank, std::size_t num, std::size_t Nx, std::size_t Ny, std::size_t steps,
    std::string const& output_name)
{

    std::vector<std::vector<double>> U(2);

    auto numa_domains = hpx::compute::host::numa_domains();

    U[0].resize(Nx * Ny, 0.0);
    U[1].resize(Nx * Ny, 0.0);

    executor_type executor(numa_domains);
    // default prefetching distance factor
    std::size_t prefetch_distance_factor = 1;
    
    auto policy = hpx::parallel::par.on(executor);

    // Iterate over the interior: skip the last and first element 
    auto f = [&](int i) { 
        line_update_prefetch(i + 1, i + Nx - 1, U[0], U[1]); 
    };
    
    hpx::util::high_resolution_timer t;    
    
    for (std::size_t t = 0; t < steps; ++t)
    {

        auto begin = U[0].begin() + Nx * t;
        auto end = U[0].begin() + Nx * t + Ny - 1;

        // Update our interior spatial domain        
        
//DETERMING PREFETCHER DISTANCE BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, hpx::parallel::prefetching_distance_determination({hpx::get_os_thread_count(), 3, 0, 0, std::size_t(std::distance(begin, end)), 0}), U[0], U[1]), 
           begin, end, f);
    }    
    double elapsed = t.elapsed();

    if (rank == 0)
        std::cout << "\n";

    if (rank == 0)
    {
        double mlups = (((Nx - 2.) * (Ny * num - 2.) * steps) / 1e6)/ elapsed;
        std::cout << "MLUPS: " << mlups << "\n";
    }
}

///////////////////////////////////////////////////////////////////////////////

int hpx_main(boost::program_options::variables_map& vm)
{
    std::size_t Nx = vm["Nx"].as<std::size_t>();
    std::size_t Ny_global = vm["Ny"].as<std::size_t>();
    std::size_t steps = vm["steps"].as<std::size_t>();

    std::size_t rank = hpx::get_locality_id();
    std::size_t num_localities = hpx::get_num_localities(hpx::launch::sync);
    std::size_t num_local_partitions = vm["local-partitions"].as<std::size_t>();

    std::size_t num_partitions = num_localities * num_local_partitions;

    // We divide our grid in stripes along the y axis.
    std::size_t Ny = Ny_global / num_partitions;

    ////////////////////////////////////////////////////////////////////////
    
    std::size_t t_origin = hpx::util::high_resolution_clock::now();

    std::vector<hpx::future<void>> workers;
    workers.reserve(num_local_partitions);
    for (std::size_t part = 0; part != num_local_partitions; ++part)
    {
        std::string output_name;
        if (vm.count("output"))
            output_name = vm["output"].as<std::string>();

        workers.push_back(hpx::async(&worker_original,
            (rank * num_local_partitions) + part, num_partitions,
            Nx, Ny, steps,
            output_name
        ));
    }

    hpx::when_all(workers).get();

    std::size_t elapsed_origin = hpx::util::high_resolution_clock::now() - t_origin;

    ////////////////////////////////////////////////////////////////////////
    // [2] Efficient chunk size
        std::size_t t_chunk = hpx::util::high_resolution_clock::now();

    for (std::size_t part = 0; part != num_local_partitions; ++part)
    {
        std::string output_name;
        if (vm.count("output"))
            output_name = vm["output"].as<std::string>();

        worker_chunk((rank * num_local_partitions) + part, 
                        num_partitions,
                        Nx, Ny, steps,
                        output_name);
    }

    std::size_t elapsed_chunk = hpx::util::high_resolution_clock::now() - t_chunk;

    ////////////////////////////////////////////////////////////////////////
    // [3] Prefetching:
    
    std::size_t t_prefetch = hpx::util::high_resolution_clock::now();

    for (std::size_t part = 0; part != num_local_partitions; ++part)
    {
        std::string output_name;
        if (vm.count("output"))
            output_name = vm["output"].as<std::string>();
        
        worker_prefecth((rank * num_local_partitions) + part, 
                        num_partitions,
                        Nx, Ny, steps,
                        output_name);
    }

    std::cout << "c\n";

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
        ("Nx", value<std::size_t>()->default_value(1024),
         "Elements in the x direction")
        ("Ny", value<std::size_t>()->default_value(1024),
         "Elements in the y direction")
        ("steps", value<std::size_t>()->default_value(100),
         "Number of steps to apply the stencil")
        ("local-partitions", value<std::size_t>()->default_value(1),
         "Number of local partitions on one locality")
        ("output", value<std::string>(),
         "Save output to file")
    ;

    // Initialize and run HPX, this example requires to run hpx_main on all
    // localities
    std::vector<std::string> const cfg = {
        "hpx.run_hpx_main!=1",
        "hpx.numa_sensitive=2"
    };

    return hpx::init(desc_commandline, argc, argv, cfg);
}
