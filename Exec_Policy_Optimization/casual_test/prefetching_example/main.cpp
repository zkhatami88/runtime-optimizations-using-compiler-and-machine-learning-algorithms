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


int hpx_main(int argc, char* argv[])
{
    // prefetching:
    std::size_t prefetch_distance_factor = 2;
    std::vector<double> c(10007, 1.0);
    std::vector<double> d(10007, 1.0);

    std::vector<std::size_t> range(10007);
    std::iota(range.begin(), range.end(), 0);

    auto policy = hpx::parallel::par;

    auto f = [&](std::size_t i) {
        c[i] = 42.1; 
    };

    // New Method
    /*
    hpx::util::high_resolution_timer t;    
    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(prefetch_distance_factor, c, d), c.begin(), c.end(), f);
    std::cout << "\n time :" << t.elapsed() << std::endl;
  
    // Old Method
    /*
    hpx::util::high_resolution_timer t;
    auto ctx = hpx::parallel::util::make_prefetcher_context(
        range.begin(), range.end(), prefetch_distance_factor, c);    
    hpx::parallel::for_each(policy, ctx.begin(), ctx.end(), f);
    std::cout << "\n time :" << t.elapsed() << std::endl;
    */

    // new exec wrapper
    hpx::parallel::parallel_executor par_exec_;
    hpx::parallel::sequential_executor seq_exec_;
    //hpx::parallel::for_each(hpx::parallel::seq.on(par_exec_), c.begin(), c.end(), f); NOT Working
    hpx::parallel::for_each(hpx::parallel::seq.on(hpx::parallel::seq_wrapper(par_exec_)), c.begin(), c.end(), f);
    hpx::parallel::for_each(hpx::parallel::seq.on(hpx::parallel::seq_wrapper(seq_exec_)), c.begin(), c.end(), f);

    return hpx::finalize();
}


int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}