#include <limits>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include "regression_models.hpp"

int main(int argc, const char * argv[]) {
	std::size_t number_of_experiments_two_class = 41;
	std::size_t number_of_features_two_class = 4;

	float** experimental_results;
	int* targets;
	float threshold = 0.15;

	

    //reading data input from .txt file
    std::string::size_type sz;
    std::string line;
    //for two classes
    //std::ifstream myfile ("training_data_two_class.txt");
    //for four classes
    //std::ifstream myfile ("training_data_multi_class.txt");
    
    
    

    /*
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
    */
	
    //learning two classes 
    /*      
    std::cout <<"\nBinary logistic regression model : \n"<<std::endl;
	learning_binary_regression_model my_nw(number_of_experiments, number_of_features_two_class,
						                  threshold, experimental_results, targets);

	my_nw.learning_two_classes();
	std::cout<<"Learning has been done!\n"<<std::endl;
    //my_nw.printing_predicted_output_two_class();
    */

    //learning multi classes
    std::cout <<"\nMulti-class logistic regression model : \n"<<std::endl;
    std::size_t number_of_experiments_multi_class = 103;
    std::size_t number_of_features_multi_class = 3;
    std::size_t number_of_multi_classes = 3;
    MatrixXf experimental_results_multi_class = MatrixXf::Random(number_of_experiments_multi_class, number_of_features_multi_class);

    targets = new int[number_of_experiments_multi_class];

    //for three classes
    std::ifstream myfile ("three_class2.txt");

    std::size_t e = 0;
    while(e < number_of_experiments_multi_class) {
        getline(myfile, line);
        std::stringstream ss(line);
        std::string str;
        std::size_t f = 0;
        while(getline(ss, str, '\t')) {
            if(f < number_of_features_multi_class) {
                experimental_results_multi_class(e, f) = std::atof(str.c_str());
            }
            else {
                targets[e] = std::stoi(str, &sz);
            }
            f++;
        }
        e++;    
    }
/*
    std::cout <<"\n ====================\n"<<std::endl;
    std::cout<<"The values of training data are : "<<std::endl;
    //computing average and variance values for each feature
    for(std::size_t i = 0; i < number_of_experiments_multi_class; i++) {        
        for(std::size_t j = 0; j < number_of_features_multi_class; j++) {
            std::cout<<experimental_results_multi_class(i, j)<<"\t";
        }
        std::cout<<targets[i]<<std::endl;
    }
    std::cout <<"\n ====================\n"<<std::endl;
     */
    multinomial_regression_model_gradient_descent my_nw(number_of_experiments_multi_class, number_of_features_multi_class,
                              number_of_multi_classes, threshold, experimental_results_multi_class, targets);

    my_nw.learning_multi_classes();
    std::cout<<"\nLearning has been done!\n"<<std::endl;
    std::cout<<"\nThe predicated class for each experimental results are: "<<std::endl;
    //my_nw.printing_predicted_output_multi_class();

	return 0;
}