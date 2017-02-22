#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>


#define num_iters 100

#define num_rows 100

int hpx_main(int argc, char* argv[])
{
    auto r = boost::irange(0, num_iters);
        
    auto m = std::vector<std::vector<double > >(num_rows);
    for (auto& v : m)
    {
        v = std::vector<double>(num_rows);
        std::generate(v.begin(), v.end(), std::rand);
    }

    auto v = std::vector<double>(num_rows);
    std::generate(v.begin(), v.end(), std::rand);
    
    std::vector<double> ret(num_rows, 0);

    auto f = [&](unsigned i)
    {        
        for (int j = 0; j < num_rows; ++j)
            ret[i] += m[i][j] * v[j];
    };

    std::cout << hpx::get_os_thread_count() << " ";
    std::cout << std::distance(r.begin(), r.end()) << " ";
   
    auto policy = hpx::parallel::par;
    std::size_t dist1 = 1;
    std::size_t dist2 = 5;
    std::size_t dist3 = 10;
    std::size_t dist4 = 50;
    std::size_t dist5 = 100;
    std::size_t dist6 = 500;
    
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist1, ret, m, v), r.begin(), r.end(), f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist2, ret, m, v), r.begin(), r.end(), f);

    std::cout << t2.elapsed() << " ";

    hpx::util::high_resolution_timer t3;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist3, ret, m, v), r.begin(), r.end(), f);

    std::cout << t3.elapsed() << " ";

    hpx::util::high_resolution_timer t4;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist4, ret, m, v), r.begin(), r.end(), f);

    std::cout << t4.elapsed() << " ";

    hpx::util::high_resolution_timer t5;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist5, ret, m, v), r.begin(), r.end(), f);

    std::cout << t5.elapsed() << " ";

    hpx::util::high_resolution_timer t6;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist6, ret, m, v), r.begin(), r.end(), f);

    std::cout << t6.elapsed() << std::endl;

    return hpx::finalize();
}