/*
#include <iostream>

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
#include <hpx/parallel/executors/adaptive_chunk_size.hpp>


int hpx_main(int argc, char* argv[])
{
    // prefetching:
    std::size_t prefetch_distance_factor = 2;
    std::vector<double> c(10007, 1.0);
    std::vector<double> d(10007, 1.0);

    std::vector<std::size_t> range(10007);
    std::iota(range.begin(), range.end(), 0);

    hpx::parallel::parallel_executor par_exec_;
    hpx::parallel::sequential_executor seq_exec_;
    auto policy = hpx::parallel::par.with(par_exec_);

    int upper_bound = 100;
    auto f = [&](std::size_t i) {
        //for(int i = 0; i < 100; i++) {
            c[i] = 42.1;
        //}         
    };

    // New Method
    
    //hpx::util::high_resolution_timer t;    
    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, prefetch_distance_factor, c, d), 
                            c.begin(), c.end(), f);

    //hpx::parallel::for_each(hpx::parallel::par.with(hpx::parallel::adaptive_chunk_size()), c.begin(), c.end(), f);
    //std::cout << "\n time :" << t.elapsed() << std::endl;
  
    // Old Method
    
    //hpx::util::high_resolution_timer t;
    //auto ctx = hpx::parallel::util::make_prefetcher_context(
    //    range.begin(), range.end(), prefetch_distance_factor, c);    
    //hpx::parallel::for_each(policy, ctx.begin(), ctx.end(), f);
    //std::cout << "\n time :" << t.elapsed() << std::endl;
    

    // new exec wrapper
    
    //hpx::parallel::for_each(hpx::parallel::seq.on(par_exec_), c.begin(), c.end(), f); NOT Working
    //hpx::parallel::for_each(hpx::parallel::seq.on(hpx::parallel::seq_wrapper(par_exec_)), c.begin(), c.end(), f);
    //hpx::parallel::for_each(hpx::parallel::seq.on(hpx::parallel::seq_wrapper(seq_exec_)), c.begin(), c.end(), f);

    // new prefetching param
    //hpx::parallel::prefetching_parameters pref_param(prefetch_distance_factor);
    //hpx::parallel::for_each(hpx::parallel::seq.with(pref_param), c.begin(), c.end(), f);
    //hpx::parallel::for_each(hpx::parallel::seq.with(prefetch_distance_factor, c), c.begin(), c.end(), f);

    return hpx::finalize();
}


int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
*/

#include <hpx/hpx_init.hpp>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>

#define size_vector 7000

int hpx_main(int argc, char* argv[])
{
    std::vector<double> a(size_vector, 1.5);
    std::vector<double> b(size_vector, 1.5);
    std::vector<double> c(size_vector, 1.5);
    
    auto f = [&](std::size_t i)
    {
        c[i] = a[i] + b[i];

        c[i] += 2.5;

        b[i] = (c[i] * a[i] + 100)/c[i];

        c[i] = a[i] * 9.9 + a[i] + 12 - b[i]/14 + b[i] * a[i];
    };

    std::cout << hpx::get_os_thread_count() << " ";
    std::cout << std::distance(a.begin(), a.end()) << " ";
   
    auto policy = hpx::parallel::par;
    std::size_t dist1 = 5;
    std::size_t dist2 = 50;
    std::size_t dist3 = 100;
    std::size_t dist4 = 500;
    std::size_t dist5 = 1000;
    std::size_t dist6 = 5000;
    
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist1, a, b, c), a.begin(), a.end(), f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist2, a, b, c), a.begin(), a.end(), f);

    std::cout << t2.elapsed() << " ";

    hpx::util::high_resolution_timer t3;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist3, a, b, c), a.begin(), a.end(), f);

    std::cout << t3.elapsed() << " ";

    hpx::util::high_resolution_timer t4;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist4, a, b, c), a.begin(), a.end(), f);

    std::cout << t4.elapsed() << " ";

    hpx::util::high_resolution_timer t5;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist5, a, b, c), a.begin(), a.end(), f);

    std::cout << t5.elapsed() << " ";

    hpx::util::high_resolution_timer t6;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist6, a, b, c), a.begin(), a.end(), f);

    std::cout << t6.elapsed() << std::endl;

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    
    // Initialize HPX, run hpx_main as the first HPX thread, and
    // wait for hpx::finalize being called.
    return hpx::init(argc, argv);
}