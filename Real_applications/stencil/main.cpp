// Copyright (c) 2016 Thomas Heller
// Copyright (c) 2017 Zahra Khatami 

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
    std::cout << "Running efficient-prefetcher-example using " << num << " Partitions\n";

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
        std::cout << "MLUPS_prefetch: " << mlups << "\n";
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
