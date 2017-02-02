#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <hpx/parallel/algorithms/for_each.hpp>

#define num_iters1_mm 1000
#define num_iters2_mm 2000
#define num_iters3_mm 10000
#define num_iters4_mm 10000
#define num_rows_mm 1000


#define num_iters1_dp 1000
#define num_iters2_dp 2000
#define num_iters3_dp 500
#define num_iters4_dp 700
#define num_rows1_dp 1000
#define num_rows2_dp 500

#define num_iters1_mv 100
#define num_iters2_mv 200
#define num_iters3_mv 500
#define num_iters4_mv 1000
#define num_rows1_mv 100
#define num_rows2_mv 200
#define num_rows3_mv 500


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

    //first test:::
    auto r1_mm = boost::irange(0, num_iters1_mm);
    auto r2_mm = boost::irange(0, num_iters2_mm);
    auto r3_mm = boost::irange(0, num_iters3_mm);
    auto r4_mm = boost::irange(0, num_iters4_mm);
        
    auto m1_mm = std::vector<std::vector<double > >(num_rows_mm);
    for (auto& v : m1_mm)
    {
        v = std::vector<double>(num_rows_mm);
        std::generate(v.begin(), v.end(), std::rand);
    }    
    
    auto m2_mm = std::vector<std::vector<double > >(num_rows_mm);
    for (auto& v : m2_mm)
    {
        v = std::vector<double>(num_rows_mm);
        std::generate(v.begin(), v.end(), std::rand);
    }   

    auto ret_mm = std::vector<std::vector<double > >(num_rows_mm);
    for (auto& v : ret_mm)
    {
        v = std::vector<double>(num_rows_mm);
    }
    
        
    auto f1_mm = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        
        for (int k = 0; k < num_rows_mm; ++k)
            for (int j = 0; j < num_rows_mm; ++j)
                ret_mm[j][k] += m1_mm[j][k] * m2_mm[j][k];

        if(ret_mm[0][0]  > 0) {
            ret_mm[0][0] += 1;
        }
        else {
            ret_mm[0][0] += 2;
        }

        return ret_mm;
    };

    auto f2_mm = [&](unsigned i) -> std::vector<std::vector<double> >
    {                
        for(int h = 0; h < 100; h++)
            for (int k = 0; k < 800 ; ++k)
                for (int j = 0; j < 800 ; ++j)
                    ret_mm[j][k] += m1_mm[j][k] * m2_mm[j][k] * 100 / ret_mm[1][1];

        if(ret_mm[0][0]  > 0) {
            ret_mm[0][0] += 1;
        }
        else if(ret_mm[0][0] == 0) {
            ret_mm[0][0] += 2;
        }
        else{
            ret_mm[0][0] +=3;
        }

        if(ret_mm[1][1]  > 0) {
            ret_mm[1][1] += 1;
        }
        else if(ret_mm[0][1] == 0) {
            ret_mm[1][1] += 2;
        }
        else{
            ret_mm[1][1] +=3;
        }        

        return ret_mm;
    };

    auto f3_mm = [](unsigned i) -> unsigned
    {
        unsigned sum = 0;
        for (int j = 0; j < 10000; ++j)
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

    auto f4_mm = [](unsigned i) -> unsigned
    {
        unsigned sum = 0;
        for(int i = 0; i < 10; i++)
            for (int j = 0; j < 10000; ++j)
                sum = sum + j * i - 200 + 0.1 * j * j;

        return sum;
    };
   
    //reading weights from txt:
    std::vector<double> w = read_weights_binary_rm();  

    hpx::util::high_resolution_timer t1_mm;
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 4002007, 2000003, 1001003, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r1_mm.begin(), r1_mm.end(), f1_mm);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 384160216, 256000010, 64080107, 3}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r2_mm.begin(), r2_mm.end(), f2_mm);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 40007, 0, 10003, 1}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r3_mm.begin(), r3_mm.end(), f3_mm);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 900024, 300000, 100012, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r4_mm.begin(), r4_mm.end(), f4_mm);

    std::cout << "\n test 1 :" << t1_mm.elapsed() << std::endl;


    //test2::
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
     
        return sum;
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
     
        return sum;
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

        return ret_dp;
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

        return ret_dp;
    };

    hpx::util::high_resolution_timer t1_dp;
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 4003, 2000, 1001, 1}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r1_dp.begin(), r1_dp.end(), f1_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 4002010, 2000005, 1001004, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r2_dp.begin(), r2_dp.end(), f2_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 1001011, 500005, 250505, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r3_dp.begin(), r3_dp.end(), f3_dp);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 100100216, 50000010, 25050107, 3}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r4_dp.begin(), r4_dp.end(), f4_dp);

    std::cout << "\n test2 : " << t1_dp.elapsed() << std::endl;


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
            for (int k = 0; k < 400; ++k)
                for (int j = 0; j < 400; ++j)
                    ret3_mv[j][k] += m3_mv[j][k] * m3_mv[j][k];

        return ret3_mv;
    };

    hpx::util::high_resolution_timer t1_mv;
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 406, 203, 103, 1}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r1_mv.begin(), r1_mv.end(), f1_mv);
    
    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 8032, 4006, 2016, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r2_mv.begin(), r2_mv.end(), f2_mv);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 1001014, 500010, 250506, 2}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r3_mv.begin(), r3_mv.end(), f3_mv);

    hpx::parallel::adaptive_for_each(hpx::parallel::features_container<double>({0, 0, 1922412, 960000, 481206, 3}), hpx::parallel::weights_container<double>({w[0], w[1], w[2], w[3], w[4], w[5], w[6]}), 
                            hpx::parallel::par, r4_mv.begin(), r4_mv.end(), f4_mv);;

    std::cout << "\n test 3 : " << t1_mv.elapsed() << std::endl;


    return hpx::finalize();
}