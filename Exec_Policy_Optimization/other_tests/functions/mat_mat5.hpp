#include <vector>
#include <algorithm>


#define num_iters 1000

#define num_rows 1000

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
        {
            for (int j = 0; j < num_rows; ++j)
                ret[i][k] += m1[i][j] * m2[j][k];
            if (ret[i][k] > 0)
                ret[i][k] * = 4.3;
        }

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