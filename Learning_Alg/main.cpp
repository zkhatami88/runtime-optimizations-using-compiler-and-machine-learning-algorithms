#include <limits>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include "regression_models.hpp"

void reading_input_values(std::size_t number_of_experiments, std::size_t number_of_features, float** experimental_results, int* targets, std::ifstream& myfile) {
    std::string::size_type sz;
    std::string line;
    std::size_t e = 0;    

    while(e < number_of_experiments) {
        getline(myfile, line);
        std::stringstream ss(line);
        std::string str;
        std::size_t f = 0;
        while(getline(ss, str, '\t')) {
            if(f < number_of_features) {
                experimental_results[e][f] = std::atof(str.c_str());
            }
            else {
                targets[e] = std::stoi(str, &sz);
            }
            f++;
        }
        e++;    
    }
}

void printing_input_data_value(std::size_t number_of_experiments, std::size_t number_of_features, float** experimental_results, int* targets) {
    std::cout <<"\n ====================\n"<<std::endl;
    std::cout<<"The values of training data are : "<<std::endl;
    for(std::size_t i = 0; i < number_of_experiments; i++) {        
        for(std::size_t j = 0; j < number_of_features; j++) {
            std::cout<<experimental_results[i][j]<<"\t";
        }
        std::cout<<targets[i]<<std::endl;
    }
    std::cout <<"\n ====================\n"<<std::endl;
}

int main(int argc, const char * argv[]) {	
	float threshold = 0.05;   

    //learning two classes    
    std::cout <<"\nBinary logistic regression model : \n"<<std::endl;

    //number of experiments and number of feautures in each experiments:
    std::size_t number_of_experiments_two_class = 41;
    std::size_t number_of_features_two_class = 4;

    //initializing
    float** experimental_results_two_class = new float*[number_of_experiments_two_class];
    for(std::size_t n = 0; n < number_of_experiments_two_class; n++) {
        experimental_results_two_class[n] = new float[number_of_features_two_class];
    }
    int* targets_two_class = new int[number_of_experiments_two_class];    

    //reading input data
    std::ifstream myfile ("training_data_two_class.txt");
    reading_input_values(number_of_experiments_two_class, number_of_features_two_class, experimental_results_two_class, targets_two_class, myfile);

    //printing inpput data values
    printing_input_data_value(number_of_experiments_two_class, number_of_features_two_class, experimental_results_two_class, targets_two_class);
    
	learning_binary_regression_model my_nw(number_of_experiments_two_class, number_of_features_two_class,
						                  threshold, experimental_results_two_class, targets_two_class);

	my_nw.learning_two_classes();
	std::cout<<"Learning has been done!\n"<<std::endl;
    std::cout<<"\nThe predicated weights for each features are: "<<std::endl;
    MatrixXf weights = my_nw.retrieving_weights_two_classes();
    my_nw.printing_predicted_output_two_class();
    
    /*
    //learning multi classes    
    std::cout <<"\nMulti-class logistic regression model : \n"<<std::endl;

    //number of experiments, number of feautures in each experiments and number of classes:
    std::size_t number_of_experiments_multi_class = 103;
    std::size_t number_of_features_multi_class = 3;
    std::size_t number_of_multi_classes = 3;

    //initializing    
    int* targets_multi_class = new int[number_of_experiments_multi_class];
    float** experimental_results_multi_class;
    experimental_results_multi_class = new float*[number_of_experiments_multi_class];
    for(std::size_t n = 0; n < number_of_experiments_multi_class; n++) {
        experimental_results_multi_class[n] = new float[number_of_features_multi_class];
    }    

    //reading input data
    std::ifstream myfile ("three_class2.txt");
    reading_input_values(number_of_experiments_multi_class, number_of_features_multi_class, experimental_results_multi_class, targets_multi_class, myfile);

    //printing inpput data values
    printing_input_data_value(number_of_experiments_multi_class, number_of_features_multi_class, experimental_results_multi_class, targets_multi_class);   


    //[1] Gradient Descent:     
    multinomial_regression_model_gradient_descent my_nw(number_of_experiments_multi_class, number_of_features_multi_class,
                                                        number_of_multi_classes, threshold, experimental_results_multi_class, targets_multi_class);    
    
    //[2] Newton-Raphson
    
    //multinomial_regression_model_newton_raphson my_nw(number_of_experiments_multi_class, number_of_features_multi_class,
    //                                                    number_of_multi_classes, threshold, experimental_results_multi_class, targets_multi_class);
    

    my_nw.learning_multi_classes();
    std::cout<<"\nLearning has been done!\n"<<std::endl;
    std::cout<<"\nThe predicated weights for each features are: "<<std::endl;
    MatrixXf weights = my_nw.retrieving_weights_multi_classes();
    my_nw.printing_predicted_output_multi_class();
    */

	return 0;
}