#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>


#define num_iters1_mv 10//100
#define num_iters2_mv 10//200
#define num_iters3_mv 10//500
#define num_iters4_mv 10//1000

#define num_rows1_mv 10//100
#define num_rows2_mv 10//200
#define num_rows3_mv 10//500

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
    auto r1_mv = boost::irange(0, num_iters1_mv);
    auto r2_mv = boost::irange(0, num_iters2_mv);
    auto r3_mv = boost::irange(0, num_iters3_mv);
    auto r4_mv = boost::irange(0, num_iters4_mv);
        
    auto m_mv = std::vector<std::vector<double > >(num_rows1_mv);
    for (auto& v : m_mv)
    {
        v = std::vector<double>(num_rows1_mv);
        std::generate(v.begin(), v.end(), std::rand);
    }

    auto v_mv = std::vector<double>(num_rows1_mv);
    std::generate(v_mv.begin(), v_mv.end(), std::rand);
    
    std::vector<double> ret1_mv_mv(num_rows1_mv, 0);

    auto f1_mv = [&](unsigned i) -> std::vector<double>
    {
        
        for (int j = 0; j < num_rows1_mv; ++j)
            ret1_mv_mv[j] += m_mv[j][j] * v_mv[j];

        if(ret1_mv_mv[0] == 0) {
            ret1_mv_mv[0]++;
        }

        if(ret1_mv_mv[2] == 0) {
            ret1_mv_mv[2] = 1;
        }

        return ret1_mv_mv;
    };

    auto m2_mv = std::vector<std::vector<double > >(num_rows2_mv);
    for (auto& v : m2_mv)
    {
        v = std::vector<double>(num_rows2_mv);
        std::generate(v.begin(), v.end(), std::rand);
    }

    auto v2_mv = std::vector<double>(num_rows2_mv);
    std::generate(v2_mv.begin(), v2_mv.end(), std::rand);
    
    std::vector<double> ret2_mv(num_rows2_mv, 0);

    auto f2_mv = [&](unsigned i) -> std::vector<double>
    {
        for(int h = 0; h < 10; h++)
            for (int j = 0; j < num_rows2_mv; ++j)
                ret2_mv[j] += m2_mv[j][j] * v2_mv[j];

        if(ret2_mv[0] == 0) {
            ret2_mv[0]++;
        }

        if(ret2_mv[2] == 0) {
            ret2_mv[2] = 1;
        }

        if(ret2_mv[3] == 0) {
            ret2_mv[3]++;
        }

        if(ret2_mv[4] == 0) {
            ret2_mv[4] = 1;
        }

        return ret2_mv;
    };

    auto m3_mv = std::vector<std::vector<double > >(num_rows3_mv);
    for (auto& v : m3_mv)
    {
        v = std::vector<double>(num_rows3_mv);
        std::generate(v.begin(), v.end(), std::rand);
    }   

    auto ret3_mv = std::vector<std::vector<double > >(num_rows3_mv);
    for (auto& v : ret3_mv)
    {
        v = std::vector<double>(num_rows3_mv);
    }
    
        
    auto f3_mv = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        
        for (int k = 0; k < num_rows3_mv; ++k)
            for (int j = 0; j < num_rows3_mv; ++j)
                ret3_mv[j][k] += m3_mv[j][k] * m3_mv[j][k];

        if(ret3_mv[0][0]  > 0) {
            ret3_mv[0][0] += 1;
        }
        else if(ret3_mv[0][0] == 0) {
            ret3_mv[0][0] += 2;
        }
        else{
            ret3_mv[0][0] +=3;
        }

        if(ret3_mv[1][1]  > 0) {
            ret3_mv[1][1] += 1;
        }
        else if(ret3_mv[0][1] == 0) {
            ret3_mv[1][1] += 2;
        }
        else{
            ret3_mv[1][1] +=3;
        }

        return ret3_mv;
    };

    auto f4_mv = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        for(int h = 0; h < 3; h++)
            for (int k = 0; k < num_rows3_mv /*- 100*/; ++k)
                for (int j = 0; j < num_rows3_mv /*- 100*/; ++j)
                    ret3_mv[j][k] += m3_mv[j][k] * m3_mv[j][k];

        return ret3_mv;
    };
   
    //reading weights from txt:
    std::vector<double> w = read_weights_binary_rm();

    hpx::util::high_resolution_timer t1_mv;
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r1_mv.begin(), r1_mv.end(), f1_mv);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r2_mv.begin(), r2_mv.end(), f2_mv);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r3_mv.begin(), r3_mv.end(), f3_mv);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::seq, r4_mv.begin(), r4_mv.end(), f4_mv);;

    std::cout << t1_mv.elapsed() << std::endl;

    return hpx::finalize();
}