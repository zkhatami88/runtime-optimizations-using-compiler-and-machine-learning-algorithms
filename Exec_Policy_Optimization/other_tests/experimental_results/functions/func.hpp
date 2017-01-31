#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>

#define num_iters1 10//1000
#define num_iters2 20//2000
#define num_iters3 10//10000000
#define num_iters4 10//10000000;

#define num_rows 10//1000

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

    return weights;
}


int hpx_main(int argc, char* argv[])
{
    auto r1 = boost::irange(0, num_iters1);
    auto r2 = boost::irange(0, num_iters2);
    auto r3 = boost::irange(0, num_iters3);
    auto r4 = boost::irange(0, num_iters4);
        
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
    
        
    auto f1 = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        
        for (int k = 0; k < num_rows; ++k)
            for (int j = 0; j < num_rows; ++j)
                ret[j][k] += m1[j][k] * m2[j][k];

        if(ret[0][0]  > 0) {
            ret[0][0] += 1;
        }
        else {
            ret[0][0] += 2;
        }

        return ret;
    };

    auto f2 = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        for(int h = 0; h < 100; h++)
            for (int k = 0; k < num_rows /*- 100*/; ++k)
                for (int j = 0; j < num_rows /*- 100*/; ++j)
                    ret[j][k] += m1[j][k] * m2[j][k] * 100 / ret[1][1];

        if(ret[0][0]  > 0) {
            ret[0][0] += 1;
        }
        else if(ret[0][0] == 0) {
            ret[0][0] += 2;
        }
        else{
            ret[0][0] +=3;
        }

        if(ret[1][1]  > 0) {
            ret[1][1] += 1;
        }
        else if(ret[0][1] == 0) {
            ret[1][1] += 2;
        }
        else{
            ret[1][1] +=3;
        }        

        return ret;
    };

    auto f3 = [](unsigned i) -> unsigned
    {
        unsigned sum = 0;
        for (int j = 0; j < 10/*100000*/; ++j)
            sum += j * i;
 
        if(sum  > 0) {
            sum += 1;
        }
        else if(sum == 0) {
            sum += 2;
        }
        else{
            sum +=3;
        }

        return sum;
    };

    auto f4 = [](unsigned i) -> unsigned
    {
        unsigned sum = 0;
        for(int i = 0; i < 10; i++)
            for (int j = 0; j < 10/*1000000*/; ++j)
                sum = sum + j * i - 200 + 0.1 * j * j;

        return sum;
    };
   
    //reading weights from txt:
    std::vector<double> w = read_weights_binary_rm();  

    hpx::util::high_resolution_timer t1;
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 427, 203, 113, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r1.begin(), r1.end(), f1);
    /*
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r2.begin(), r2.end(), f2);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r3.begin(), r3.end(), f3);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r4.begin(), r4.end(), f4);*/

    std::cout << t1.elapsed() << std::endl;

    return hpx::finalize();
}