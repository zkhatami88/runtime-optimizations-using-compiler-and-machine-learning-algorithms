#include <limits>
#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>
#include "regression_models.hpp"

void reading_random_input_values(std::size_t number_of_experiments, std::size_t number_of_features, float** experimental_results, int* targets, std::ifstream& myfile) {    
    std::string line;
    std::size_t e = 0;

    //first line contains only name of features for each experiment:   
    getline(myfile, line);

    //assigning values for experimental_results:
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
                targets[e] = std::stoi(str);
            }
            f++;
        }
        e++;    
    }
}

void reading_real_input_values(std::size_t number_of_experiments, std::size_t number_of_features, float** experimental_results, int* targets, std::ifstream& myfile) {
    std::string line;
    std::size_t e = 0;
    std::size_t num_seq = 0;
    std::size_t num_par = 0;

    //first line contains only name of features for each experiment:   
    getline(myfile, line);

    //assigning values for experimental_results:
    while(e < number_of_experiments) {
        getline(myfile, line);
        std::stringstream ss(line);
        std::string str;
        std::size_t f = 0;

        //first feature is file_name
        getline(ss, str, ' ');

        //reading features:
        while(f < number_of_features) {
            getline(ss, str, ' ');
            experimental_results[e][f] = std::atof(str.c_str());
            f++;
        }

        //skipping last three features:
        getline(ss, str, ' ');
        getline(ss, str, ' ');
        getline(ss, str, ' ');
        getline(ss, str, ' '); //:nan : as these feature is the same as prev one
        //getline(ss, str, ' '); // : 9%
        //getline(ss, str, ' '); //: 30%

        //reading exec time for assiging target value (0 or 1)
        getline(ss, str, ' ');
        float time_seq = std::atof(str.c_str());
        getline(ss, str, ' ');
        float time_par = std::atof(str.c_str());
        if(time_seq < time_par) {
            targets[e] = 0;
            num_seq++;
        }
        else {
            targets[e] = 1;
            num_par++;
        }
        e++;    
    }

    //printing diversity of input data:
    std::cout<<"num_seq = "<<num_seq<<", num_par = "<<num_par<<std::endl;
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
	float threshold = 0.2;
    std::string line;
    
    //learning two classes    
    std::cout <<"\nBinary logistic regression model : \n"<<std::endl;

    //reading input data : number of experiments and number of feautures in each experiments
    //std::ifstream myfile ("inputs/training_data_two_class.txt");
    std::ifstream myfile ("inputs/data.dat");
    getline(myfile, line);
    std::stringstream ss(line);
    std::string str;
    getline(ss, str, ' ');
    std::size_t number_of_experiments_two_class = std::stoi(str);
    getline(ss, str, ' ');
    std::size_t number_of_features_two_class = std::stoi(str);
    
    //initializing
    float** experimental_results_two_class = new float*[number_of_experiments_two_class];
    for(std::size_t n = 0; n < number_of_experiments_two_class; n++) {
        experimental_results_two_class[n] = new float[number_of_features_two_class];
    }
    int* targets_two_class = new int[number_of_experiments_two_class];    

    //reading input data
    //reading_random_input_values(number_of_experiments_two_class, number_of_features_two_class, experimental_results_two_class, targets_two_class, myfile);

    reading_real_input_values(number_of_experiments_two_class, number_of_features_two_class, experimental_results_two_class, targets_two_class, myfile);

    //printing inpput data values
    //printing_input_data_value(number_of_experiments_two_class, number_of_features_two_class, experimental_results_two_class, targets_two_class);
    
    
	learning_binary_regression_model my_nw(number_of_experiments_two_class, number_of_features_two_class,
						                  threshold, experimental_results_two_class, targets_two_class);

	my_nw.learning_two_classes();
	std::cout<<"Learning has been done!\n"<<std::endl;
    std::cout<<"\nThe predicated weights for each features are: "<<std::endl;
    my_nw.retrieving_weights_two_classes_into_txt_file();
    my_nw.printing_predicted_output_two_class();
    
    /*
    //learning multi classes    
    std::cout <<"\nMulti-class logistic regression model : \n"<<std::endl;

    //reading input data : number of experiments and number of feautures in each experiments
    std::ifstream myfile ("inputs/three_class2.txt");
    getline(myfile, line);
    std::stringstream ss(line);
    std::string str;
    getline(ss, str, ' ');
    std::size_t number_of_experiments_multi_class = std::stoi(str);
    getline(ss, str, ' ');
    std::size_t number_of_features_multi_class = std::stoi(str);

    //number of classes:
    std::size_t number_of_multi_classes = 3;

    //initializing    
    int* targets_multi_class = new int[number_of_experiments_multi_class];
    float** experimental_results_multi_class;
    experimental_results_multi_class = new float*[number_of_experiments_multi_class];
    for(std::size_t n = 0; n < number_of_experiments_multi_class; n++) {
        experimental_results_multi_class[n] = new float[number_of_features_multi_class];
    }    

    //reading input data    
    reading_random_input_values(number_of_experiments_multi_class, number_of_features_multi_class, experimental_results_multi_class, targets_multi_class, myfile);

    //printing inpput data values
    printing_input_data_value(number_of_experiments_multi_class, number_of_features_multi_class, experimental_results_multi_class, targets_multi_class);   

    //choose one of them: Gradient Descen OR Newton-Raphson
    //[1] Gradient Descent:     
    multinomial_regression_model_gradient_descent my_nw(number_of_experiments_multi_class, number_of_features_multi_class,
                                                        number_of_multi_classes, threshold, experimental_results_multi_class, targets_multi_class);    
    
    //[2] Newton-Raphson:    
    //multinomial_regression_model_newton_raphson my_nw(number_of_experiments_multi_class, number_of_features_multi_class,
    //                                                    number_of_multi_classes, threshold, experimental_results_multi_class, targets_multi_class);    

    my_nw.learning_multi_classes();
    std::cout<<"\nLearning has been done!\n"<<std::endl;
    std::cout<<"\nThe predicated weights for each features are: "<<std::endl;
    my_nw.retrieving_weights_multi_classes_into_text_file();
    my_nw.printing_predicted_output_multi_class();
    */
	return 0;
}