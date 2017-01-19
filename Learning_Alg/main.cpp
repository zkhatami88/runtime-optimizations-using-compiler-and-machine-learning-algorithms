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
    std::size_t number_of_two_classes = 2;
    std::size_t number_of_multi_classes = 4;
	double** experimental_results;
	unsigned* targets;
	double threshold = 0.1;

	experimental_results = new double*[number_of_experiments];
	targets = new unsigned[number_of_experiments];

	for(std::size_t i = 0; i < number_of_experiments; i++) {
    	experimental_results[i] = new double[number_of_features];
    }

    //reading data input from .txt file
    std::string::size_type sz;
    std::string line;
    //for two classes
    //std::ifstream myfile ("training_data_two_class.txt");
    std::ifstream myfile ("training_data_multi_class.txt");
    
    std::size_t e = 0;
    while(e < number_of_experiments) {
    	getline(myfile, line);
    	std::stringstream ss(line);
    	std::string str;
    	std::size_t f = 0;
    	while(getline(ss, str, '\t')) {
    		if(f < number_of_features) {
    			experimental_results[e][f] = std::stod(str, &sz);
    		}
    		else {
    			targets[e] = std::stoi(str, &sz);
    		}
    		f++;
    	}
    	e++;    
    }


    std::cout <<"\n ====================\n"<<std::endl;
    std::cout<<"The values of training data are : "<<std::endl;
	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features; j++) {
			std::cout<<experimental_results[i][j]<<"\t";
		}
		std::cout<<targets[i]<<std::endl;
	}
    std::cout <<"\n ====================\n"<<std::endl;
	
    //learning two classes
    /*
    std::cout <<"\nBinary logistic regression model : \n"<<std::endl;
	learning_network my_nw(number_of_experiments, number_of_features, 
						      number_of_two_classes, threshold);

    my_nw.initilializer_two_class(experimental_results, targets);
	my_nw.learning_two_classes();
	std::cout<<"Learning has been done!\n"<<std::endl;

	double* weightsb = my_nw.retrieving_weights_two_classes();
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		std::cout<<"weights[%i] = "<<weightsb[i]<<i<<std::endl;
	}*/
    

    //learning multi classes
    //std::cout <<"\nMulti-class logistic regression model : \n"<<std::endl;
    learning_network my_nw(number_of_experiments, number_of_features, 
                              number_of_multi_classes, threshold);

    my_nw.initilializer_multi_class(experimental_results, targets);
    my_nw.learning_multi_classes();
    //std::cout<<"Learning has been done!\n"<<std::endl;

    double** weightsm = my_nw.retrieving_weights_multi_classes();
    /*
    for(std::size_t k = 0; k < number_of_multi_classes; k++) {
        for(std::size_t f = 0; f < number_of_features; f++) {
            std::size_t offset = k * number_of_features + f;
            std::cout<<"weights["<<k<<"]["<<f<<"] = "<<weightsm[offset][0]<<"\t";
        }
        std::cout<<std::endl;
    }*/

	return 0;
}