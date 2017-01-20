#include <vector>
#include <algorithm>


#define num_iters 1000

#define num_rows 1000

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

    auto f = [&](unsigned i) -> std::vector<double>
    {
        
        for (int j = 0; j < num_rows; ++j)
            ret[i] += m[i][j] * v[j];

        return ret;
    };
   
    hpx::util::high_resolution_timer t1;
    
    hpx::parallel::for_each_n(hpx::parallel::seq, r.begin(), num_iters, f);
    
    std::cout << t1.elapsed() << " ";
    
    hpx::util::high_resolution_timer t2;

    hpx::parallel::for_each_n(hpx::parallel::par, r.begin(), num_iters, f);

    std::cout << t2.elapsed() << std::endl;

    return hpx::finalize();
}