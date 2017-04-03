#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>

#define size_vector 100000

int hpx_main(int argc, char* argv[])
{
    std::vector<double> a(size_vector, 1.5);
    std::vector<double> b(size_vector, 1.5);
    std::vector<double> c(size_vector, 1.5);
    
    auto f = [&](std::size_t j)
    {
        for(int i = 0; i < 100; i++) {

            c[i] = a[i] + b[i];

            c[i] += 2.5;

            b[i] = (c[i] * a[i] + 100)/c[i];

        };
    };

    std::cout << hpx::get_os_thread_count() << " ";
    std::cout << std::distance(a.begin(), a.end()) << " ";
   
    hpx::parallel::dynamic_chunk_size dcs1(size_vector * 0.001);
    hpx::parallel::dynamic_chunk_size dcs2(size_vector * 0.01);
    hpx::parallel::dynamic_chunk_size dcs3(size_vector * 0.1);
    hpx::parallel::dynamic_chunk_size dcs4(size_vector * 1);
    hpx::parallel::dynamic_chunk_size dcs5(size_vector * 1.1);

    auto p1 = hpx::parallel::par.with(dcs1);
    auto p2 = hpx::parallel::par.with(dcs2);
    auto p3 = hpx::parallel::par.with(dcs3);
    auto p4 = hpx::parallel::par.with(dcs4);
    auto p5 = hpx::parallel::par.with(dcs5);
    
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each(p1, a.begin(), a.end(), f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each(p2, a.begin(),a.end(), f);

    std::cout << t2.elapsed() << " ";

    hpx::util::high_resolution_timer t3;

    hpx::parallel::for_each(p3, a.begin(), a.end(), f);

    std::cout << t3.elapsed() << " ";

    hpx::util::high_resolution_timer t4;

    hpx::parallel::for_each(p4, a.begin(), a.end(), f);

    std::cout << t4.elapsed() << " ";

    hpx::util::high_resolution_timer t5;

    hpx::parallel::for_each(p5, a.begin(), a.end(), f);

    std::cout << t5.elapsed() << std::endl;

    return hpx::finalize();
}