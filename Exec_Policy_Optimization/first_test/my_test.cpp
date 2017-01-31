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

#define num_iters 10

#define num_rows 10

std::vector<double> read_weights_binary_rm(std::ifstream& myfile) {
    std::vector<double> weights;
    std::string line;
    std::string str;
    //std::ifstream myfile ("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Learning_Alg/weights.txt");

    //first line contains weights:   
    getline(myfile, line);
    std::stringstream ss(line);
    while(getline(ss, str, ' ')) {
        weights.push_back(std::atof(str.c_str()));
    }

    return weights;
}

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
            for (int j = 0; j < num_rows; ++j)
                ret[i][k] += m1[i][j] * m2[j][k];

        return ret;
    };
   
   
    hpx::parallel::for_each_n(hpx::parallel::seq, r.begin(), num_iters, f);
    hpx::parallel::for_each(hpx::parallel::seq, r.begin(), r.end(), f);
    
    std::ifstream myfile ("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Learning_Alg/weights.txt");
    std::vector<double> w = read_weights_binary_rm(myfile);
    for(int i = 0; i < w.size(); i++) {
        std::cout<<w[i]<<" ";
    }
    std::cout<<std::endl;

    //these values will be overwritten by compiler
    double f0, f1, f2, f3, f4, f5;
    f0 = 7; f1 = 8; f2 = 9; f3 = 4; f4 = 5; f5 = 6;

    //initially clang assignes the value of w0 and w1 to be 0
    //which are number of threads and number of iterations
    //these values will be determined at runtime
    /*
    hpx::parallel::for_each(hpx::parallel::info_container<
                                double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6], 
    										f0, f1, f2, f3, f4, f5}), 
    						    hpx::parallel::seq, r.begin(), r.end(), f);

	/*
    hpx::parallel::info_container<int> my_info {1, 2, 3, 4, 5, 6, 7};
    auto weights = my_info.get_weights();
    auto features = my_info.get_features();

    std::cout<<"\n size of weights = "<<weights.size()<<std::endl;
    std::cout<<"\n size of features = "<<features.size()<<std::endl;

    for(auto i : weights) {
    	std::cout<<i<<", ";
    }
    std::cout<<"\n";

    for(auto i : features) {
    	std::cout<<i<<", ";
    }
    std::cout<<"\n";*/

    return hpx::finalize();
}


int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}