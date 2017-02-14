#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>

//new classes for implementing machine learning techniques
#include <hpx/parallel/seq_or_par.hpp>
#include <hpx/parallel/chunk_size_determination.hpp>
#include <hpx/parallel/prefetching_distance_determination.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>
#include <hpx/parallel/executors/adaptive_chunk_size.hpp>

#define num_iters 10

#define num_rows 10

int hpx_main(int argc, char* argv[])
{
    auto r = boost::irange(0, num_iters);
        
    auto m1 = std::vector<std::vector<double > >(num_rows);
    for (auto& v : m1)
    {
        v = std::vector<double>(num_rows);
        std::generate(v.begin(), v.end(), std::rand);
    }  
    
    auto m2 = std::vector<std::vector<double > >(num_rows);
    for (auto& v : m2)
    {
        v = std::vector<double>(num_rows);
        std::generate(v.begin(), v.end(), std::rand);
    }   

    auto ret = std::vector<std::vector<double > >(num_rows);
    for (auto& v : ret)
    {
        v = std::vector<double>(num_rows);
    }
        
    auto f = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        for (int k = 0; k < num_rows; ++k)
            for (int j = 0; j < num_rows; ++j)
                ret[i][k] += m1[i][j] * m2[j][k];

        return ret;
    };
    
    // prefetching instruction
    std::vector<double> c(10007, 1.0);
    std::vector<double> d(10007, 1.0);
    std::size_t prefetch_distance_factor = 2;
    std::vector<std::size_t> range(10007);
    std::iota(range.begin(), range.end(), 0);

    auto f_p = [&](std::size_t i) { c[i] = 42.1; };

    hpx::parallel::sequential_executor seq_exec;
    hpx::parallel::parallel_executor par_exec;
    hpx::parallel::dynamic_chunk_size dcs(10);

    //seq or par
    //DETERMING EXECUTION POLICY BASED ON STATIC AND DYNAMIC FEATURES:
 	if (hpx::parallel::seq_or_par({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})) 
 	 	hpx::parallel::for_each(hpx::parallel::seq, r.begin(), r.end(), f);
 	else 
 	 	hpx::parallel::for_each(hpx::parallel::par, r.begin(), r.end(), f);
    //DETERMING EXECUTION POLICY BASED ON STATIC AND DYNAMIC FEATURES:
 	if (hpx::parallel::seq_or_par({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})) 
 	 	hpx::parallel::for_each(hpx::parallel::seq.on(seq_exec), r.begin(), r.end(), f);
 	else 
 	 	hpx::parallel::for_each(hpx::parallel::par.on(seq_exec), r.begin(), r.end(), f);
    //DETERMING EXECUTION POLICY BASED ON STATIC AND DYNAMIC FEATURES:
 	if (hpx::parallel::seq_or_par({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})) 
 	 	hpx::parallel::for_each(hpx::parallel::seq.with(dcs), r.begin(), r.end(), f);
 	else 
 	 	hpx::parallel::for_each(hpx::parallel::par.with(dcs), r.begin(), r.end(), f);
    //DETERMING EXECUTION POLICY BASED ON STATIC AND DYNAMIC FEATURES:
 	if (hpx::parallel::seq_or_par({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})) 
 	 	hpx::parallel::for_each(hpx::parallel::seq.on(seq_exec).with(dcs), r.begin(), r.end(), f);
 	else 
 	 	hpx::parallel::for_each(hpx::parallel::par.on(seq_exec).with(dcs), r.begin(), r.end(), f);

    // choosing efficient chunk_size
    //DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::chunk_size_determination({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})),  r.begin(), r.end(), f);
    //DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::chunk_size_determination({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})).on(par_exec),  r.begin(), r.end(), f);
    //DETERMING CHUNK SIZES BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::chunk_size_determination({hpx::get_os_thread_count(), 424, 200, 112, std::size_t(std::distance(r.begin(), r.end())), 2})).on(par_exec),  r.begin(), r.end(), f);
    
    //choosing efficient prefetcher_distance_factor
    //DETERMING PREFETCHER DISTANCE BASED ON STATIC AND DYNAMIC FEATURES:
	hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(hpx::parallel::prefetching_distance_determination({hpx::get_os_thread_count(), 1, 1, 0, std::size_t(std::distance(range.begin(), range.end())), 0}), c, d),  range.begin(), range.end(), f_p);
    

    return hpx::finalize();
}
