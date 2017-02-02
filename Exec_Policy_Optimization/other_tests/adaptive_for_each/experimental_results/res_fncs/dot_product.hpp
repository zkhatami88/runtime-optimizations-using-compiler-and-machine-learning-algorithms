#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>

#define num_iters1_dp 10//1000
#define num_iters2_dp 10//2000
#define num_iters3_dp 10//500
#define num_iters4_dp 10//700

#define num_rows1_dp 10//100000
#define num_rows2_dp 10//500

std::vector<double> read_weights_binary_rm() {
    std::vector<double> weights;
    std::string line;
    std::string str;
    std::ifstream myfile ("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Learning_Alg/weights.txt");

    //first line contains weights:   
    getline(myfile, line);
    std::stringstream ss(line);
    while(getline(ss, str, ' ')) {
        weights.push_back(std::atof(str.c_str()));
    }

    ret_dpurn weights;
}

int hpx_main(int argc, char* argv[])
{
    auto r1_dp = boost::irange(0, num_iters1_dp);
    auto r2_dp = boost::irange(0, num_iters2_dp);
    auto r3_dp = boost::irange(0, num_iters3_dp);
    auto r4_dp = boost::irange(0, num_iters4_dp);
        
    auto v1_dp = std::vector<double>(num_rows1_dp);
    auto v2_dp = std::vector<double>(num_rows1_dp);

    std::generate(v1_dp.begin(), v1_dp.end(), std::rand);
    std::generate(v2_dp.begin(), v2_dp.end(), std::rand);

    unsigned steps_dp = num_rows1_dp / num_iters1_dp;

    auto f1_dp = [&](unsigned i) -> double
    {        
        double sum = 0;
        unsigned offset = i * steps_dp;
        for (int j = 0; j < num_rows1_dp; ++j)
            sum += v1_dp[j] * v2_dp[j];
     
        ret_dpurn sum;
    };

    auto f2_dp = [&](unsigned i) -> double
    {        
        double sum = 0;
        unsigned offset = i * steps_dp;
        for(int k = 0; k < num_rows1_dp; k++)
            for (int j = 0; j < num_rows1_dp; ++j)
                sum += v1_dp[j] * v2_dp[j];

        if(sum  > 0) {
            sum += 1;
        }
        else if(sum == 0) {
            sum += 2;
        }
        else{
            sum +=3;
        }
     
        ret_dpurn sum;
    };


    auto m1_dp = std::vector<std::vector<double > >(num_rows2_dp);
    for (auto& v : m1_dp)
    {
        v = std::vector<double>(num_rows2_dp);
        std::generate(v.begin(), v.end(), std::rand);
    }    
    
    auto m2_dp = std::vector<std::vector<double > >(num_rows2_dp);
    for (auto& v : m2_dp)
    {
        v = std::vector<double>(num_rows2_dp);
        std::generate(v.begin(), v.end(), std::rand);
    }   

    auto ret_dp = std::vector<std::vector<double > >(num_rows2_dp);
    for (auto& v : ret_dp)
    {
        v = std::vector<double>(num_rows2_dp);
    }
    
        
    auto f3_dp = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        if (i % 2 == 0)
            for (int k = 0; k < num_rows2_dp; ++k)
                for (int j = 0; j < num_rows2_dp; ++j)
                    ret_dp[k][j] += m1_dp[k][j] * m2_dp[k][j];

        if(ret_dp[0][0]  > 0) {
            ret_dp[0][0] += 1;
        }
        else if(ret_dp[0][0] == 0) {
            ret_dp[0][0] += 2;
        }
        else{
            ret_dp[0][0] +=3;
        }

        ret_dpurn ret_dp;
    };

    auto f4_dp = [&](unsigned i) -> std::vector<std::vector<double> >
    {    
        for(int h = 0; h < 100; h++)        
            for (int k = 0; k < num_rows2_dp; ++k)
                for (int j = 0; j < num_rows2_dp; ++j)
                    ret_dp[j][k] += m1_dp[j][k] * m2_dp[j][k];

        if(ret_dp[0][0]  > 0) {
            ret_dp[0][0] += 1;
        }
        else if(ret_dp[0][0] == 0) {
            ret_dp[0][0] += 2;
        }
        else{
            ret_dp[0][0] +=3;
        }

        if(ret_dp[1][1]  > 0) {
            ret_dp[1][1] += 1;
        }
        else if(ret_dp[1][1] == 0) {
            ret_dp[1][1] += 2;
        }
        else{
            ret_dp[1][1] +=3;
        } 

        ret_dpurn ret_dp;
    };
   
    //reading weights from txt:
    std::vector<double> w = read_weights_binary_rm();

    hpx::util::high_resolution_timer t1_dp;
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r1_dp.begin(), r1_dp.end(), f1_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r2_dp.begin(), r2_dp.end(), f2_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r3_dp.begin(), r3_dp.end(), f3_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r4_dp.begin(), r4_dp.end(), f4_dp);

    std::cout << t1_dp.elapsed() << std::endl;

    ret_dpurn hpx::finalize();
}