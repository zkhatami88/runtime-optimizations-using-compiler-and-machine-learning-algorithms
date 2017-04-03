#include <iostream>
#include <vector>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/executors/dynamic_chunk_size.hpp>


int num_iters = 900;

int num_rows = 90000;

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
        for(int h = 0; h < 2; h++)
            for(int k = 0; k < 2; k++)
                for (int j = 0; j < 90; ++j){
                    sum += v1[j] * v2[j]; 
                    if(sum == 100) 
                        sum += 50;    
                }
    };

    std::cout << hpx::get_os_thread_count() << " ";
    std::cout << std::distance(r.begin(), r.end()) << " ";

    hpx::parallel::dynamic_chunk_size dcs1(num_iters * 0.001);
    hpx::parallel::dynamic_chunk_size dcs2(num_iters * 0.01);
    hpx::parallel::dynamic_chunk_size dcs3(num_iters * 0.1);
    hpx::parallel::dynamic_chunk_size dcs4(num_iters * 1);
    hpx::parallel::dynamic_chunk_size dcs5(num_iters * 1.1);

    auto p1 = hpx::parallel::par.with(dcs1);
    auto p2 = hpx::parallel::par.with(dcs2);
    auto p3 = hpx::parallel::par.with(dcs3);
    auto p4 = hpx::parallel::par.with(dcs4);
    auto p5 = hpx::parallel::par.with(dcs5);
    
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each(p1, r.begin(), r.end(), f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each(p2, r.begin(), r.end(), f);

    std::cout << t2.elapsed() << " ";

    hpx::util::high_resolution_timer t3;

    hpx::parallel::for_each(p3, r.begin(), r.end(), f);

    std::cout << t3.elapsed() << " ";

    hpx::util::high_resolution_timer t4;

    hpx::parallel::for_each(p4, r.begin(), r.end(), f);

    std::cout << t4.elapsed() << " ";

    hpx::util::high_resolution_timer t5;

    hpx::parallel::for_each(p5, r.begin(), r.end(), f);

    std::cout << t5.elapsed() << std::endl;

    return hpx::finalize();
}