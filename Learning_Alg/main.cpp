#include <limits>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include "learning_alg.hpp"

int main(int argc, const char * argv[]) {
	std::size_t number_of_experiments = 41;
	std::size_t number_of_features = 4;
	double** experimental_results;
	double* outputs;
	double threshold = 0.1;

	experimental_results = new double*[number_of_experiments];
	outputs = new double[number_of_experiments];

	for(std::size_t i = 0; i < number_of_experiments; i++) {
    	experimental_results[i] = new double[number_of_features + 1];
    }

    //reading data input from .txt file
    std::string::size_type sz;
    std::string line;
    std::ifstream myfile ("/home/zahra/Desktop/training_data.txt");
    
    std::size_t e = 0;
    while(e < number_of_experiments) {
    	getline(myfile, line);
    	std::stringstream ss(line);
    	std::string str;
    	experimental_results[e][0] = 1;
    	std::size_t f = 1;
    	while(getline(ss, str, '\t')) {
    		if(f < number_of_features + 1) {
    			experimental_results[e][f] = std::stod(str, &sz);
    			std::cout<<experimental_results[e][f]<<"	";
    		}
    		else {
    			outputs[e] = std::stod(str, &sz);
    			std::cout<<outputs[e]<<"	";
    		}
    		f++;
    	}
    	std::cout<<std::endl;
    	e++;    
    }

    std::cout <<"\n ====================\n"<<std::endl;
    std::cout<<"The values of training data are : "<<std::endl;
	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			std::cout<<experimental_results[i][j]<<" , ";
		}
		outputs[e];
		std::cout<<std::endl;
	}
	std::cout <<"\n ====================\n"<<std::endl;
	
	learning_network my_nw(number_of_experiments, number_of_features, 
						experimental_results, outputs, threshold);


	my_nw.learning();
	std::cout <<"\n ====================\n"<<std::endl;

	std::cout<<"Learning has been done!"<<std::endl;
	double* weights = my_nw.retrieving_weights();

	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		std::cout<<"weights[%i] = "<<weights[i]<<i<<std::endl;
	}

	return 0;
}