//This file revise the generated data from creating+data test functions in a good shape

#include <math.h>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>


int main() {
	std::string line;

	// for chunk:
	//std::ofstream outputFile("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/revising_data/chunk_training_data.dat");
	//std::ifstream myfile ("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Learning_Alg/algorithms/inputs/data_chunk_3.dat");

	// for prefetch:
	std::ofstream outputFile("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Exec_Policy_Optimization/creating_data/revising_data/prefetching_training_data.dat");
	std::ifstream myfile ("/home/zahra/Desktop/runtime_opt_with_compiler_and_ML/Learning_Alg/algorithms/inputs/data_prefetch_3.dat");

	// number of expr with number of features & classes
	getline(myfile, line);
    std::stringstream ss(line);
    std::string str;
    getline(ss, str, ' ');
    std::size_t number_of_experiments_multi_class = std::stoi(str);
    getline(ss, str, ' ');
    std::size_t number_of_features_multi_class = std::stoi(str);
    getline(ss, str, ' ');
    std::size_t number_of_multi_classes = std::stoi(str);

    // printing above infor in output file
    outputFile << number_of_experiments_multi_class << " " << number_of_features_multi_class << " " << number_of_multi_classes << "\n";

    //skipping first line
    getline(myfile, line);

    //reading all experiments
    for(int i = 0; i < number_of_experiments_multi_class; i++) {
    	getline(myfile, line);
        std::stringstream ss(line);
        std::string str;

        //skipping name
        getline(ss, str, ' ');

        //skipping static informations for the first 4 loops for chunk and 5 loops for prefetch 
        for(int l = 0; l < 5; l++) {
            for(int c = 0; c < 4; c++) { //this 4 for both chunk and prefetch
                getline(ss, str, ' ');
            }
        }

        //printing features:
        int f = 0;
        while(f < number_of_features_multi_class) {
            getline(ss, str, ' ');
            outputFile << str << " ";
            f++;
        }

        //printing execution times
        for(int c = 0; c < number_of_multi_classes - 1; c++) {
            getline(ss, str, ' ');
            outputFile << str << " ";
        }

        //printing last exec time
        getline(ss, str, ' ');
        outputFile << str << "\n";
    }
	

	return 0;
}