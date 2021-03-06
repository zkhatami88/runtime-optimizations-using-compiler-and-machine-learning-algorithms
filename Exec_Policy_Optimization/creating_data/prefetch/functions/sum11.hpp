#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>

#define size_vector 10000

int hpx_main(int argc, char* argv[])
{
    auto a = std::vector<double>(size_vector);
    std::generate(a.begin(), a.end(), std::rand);

    auto b = std::vector<double>(size_vector);
    std::generate(b.begin(), b.end(), std::rand);

    auto c = std::vector<double>(size_vector);
    std::generate(c.begin(), c.end(), std::rand);
    
    auto f = [&](unsigned i)
    {
        for(int i = 0; i < size_vector; i++) {

            c[i] = a[i] + b[i];

            c[i] += 2.5;

            b[i] = (c[i] * a[i] + 100)/c[i];

            a[i] = c[i]/10 + b[i] * 8.8;

        };
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