#include <iostream>
#include <vector>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>


int num_iters = 1200;

int num_rows = 120000;

int hpx_main(int argc, char* argv[])
{
    auto r = boost::irange(0, num_iters);
        
    auto v1 = std::vector<double>(num_rows);
    auto v2 = std::vector<double>(num_rows);

    std::generate(v1.begin(), v1.end(), std::rand);
    std::generate(v2.begin(), v2.end(), std::rand);

    unsigned steps = num_rows / num_iters;

    auto f = [&](unsigned i)
    {        
        double sum = 0;
        unsigned offset = i * steps;
        for (int j = 0; j < 120; ++j)
            sum += v1[j] * v2[j]; 

        for (int j = 0; j < 120; ++j)
            sum += v1[j] * v2[j]; 

        for (int j = 0; j < 120; ++j)
            sum += v1[j] * v2[j];     
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
    
    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist1, v1, v2), r.begin(), r.end(), f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist2, v1, v2), r.begin(), r.end(), f);

    std::cout << t2.elapsed() << " ";

    hpx::util::high_resolution_timer t3;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist3, v1, v2), r.begin(), r.end(), f);

    std::cout << t3.elapsed() << " ";

    hpx::util::high_resolution_timer t4;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist4, v1, v2), r.begin(), r.end(), f);

    std::cout << t4.elapsed() << " ";

    hpx::util::high_resolution_timer t5;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist5, v1, v2), r.begin(), r.end(), f);

    std::cout << t5.elapsed() << " ";

    hpx::util::high_resolution_timer t6;

    hpx::parallel::for_each(hpx::parallel::execution::make_prefetcher_policy(policy, dist6, v1, v2), r.begin(), r.end(), f);

    std::cout << t6.elapsed() << std::endl;

    return hpx::finalize();
}