#include <vector>
#include <algorithm>

int num_iters = 1000;

int num_rows = 100000;

int hpx_main(int argc, char* argv[])
{
    auto r = boost::irange(0, num_iters);
        
    auto v1 = std::vector<double>(num_rows);
    auto v2 = std::vector<double>(num_rows);

    std::generate(v1.begin(), v1.end(), std::rand);
    std::generate(v2.begin(), v2.end(), std::rand);

    unsigned steps = num_rows / num_iters;

    auto f = [&](unsigned i) -> double
    {        
        double sum = 0;
        unsigned offset = i * steps;
        for (int j = 0; j < 100; ++j)
            sum += v1[j] * v2[j];
     
        return sum;
    };
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each_n(hpx::parallel::seq, r.begin(), num_iters, f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each_n(hpx::parallel::par, r.begin(), num_iters, f);

    std::cout << t2.elapsed() << std::endl;

    return hpx::finalize();
}