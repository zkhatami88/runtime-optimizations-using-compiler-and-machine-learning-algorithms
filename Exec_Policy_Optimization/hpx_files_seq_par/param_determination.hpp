//  Copyright (c) 2017 Zahra Khatami
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/// \file parallel/param_determination.hpp
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>

#if !defined(HPX_PARALLEL_PARAM_DETERMINATION_FEB_1_2017_0300PM)
#define HPX_PARALLEL_PARAM_DETERMINATION_FEB_1_2017_0300PM

namespace hpx { namespace parallel {
    
std::size_t param_determination(std::vector<std::size_t> features)
{   
    static std::vector<std::vector<double>> weights;
        
    if (weights.size() == 0)
    {    	

    	std::cout<<"\nPlease include normalization parameters (variance and average) in the first line \n"; 
        std::ifstream infile("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/other_tests/seq_par/test2/weights_param_determination.dat");

        // first line include (var, average) of each features
        std::string line_normalization;
        getline(infile, line_normalization);
        std::vector<double> normalization_params;
        std::istringstream normalization_buffer(line_normalization);
        std::transform(std::istream_iterator<std::string>(normalization_buffer), 
                    std::istream_iterator<std::string>(),
                    std::back_inserter(normalization_params),
                    [](std::string s) {return std::stod(s);});

        // the rest of the lines include the values of weights
        for(std::string line; getline( infile, line ); )
        {
        	std::vector<double> weights_row;
            std::istringstream buffer(line);
            std::size_t t = 0;

            std::transform(std::istream_iterator<std::string>(buffer), 
                      std::istream_iterator<std::string>(),
                      std::back_inserter(weights_row),
                      [&t, &normalization_params](std::string s) {
                      	double w = double((std::stod(s) - normalization_params[t * 2])/normalization_params[2 * t + 1]);
                      	t++;
                      	return w;
                      });
            weights.push_back(weights_row);
        }      
    }

    features[0] = hpx::get_os_thread_count();
    
    //initial class = 0
    std::size_t determined_class = 0;
    double max = 0.0;

    //max of (wi * f)
    for(std::size_t w = 0; w < weights.size(); w++) {
    	double sum = 0.0;
    	for(std::size_t f = 0; f < features.size(); f++) {
    		sum += weights[w][f] * features[f];
    	}
    	if(max < sum) {
    		max = sum;
    		determined_class = w;
    	}
    }

    return determined_class;
}

}}

#endif
