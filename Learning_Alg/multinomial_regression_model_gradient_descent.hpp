#include <limits>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/LU>

using namespace Eigen;

#define MAX_FLOAT (std::numeric_limits<float>::max())
#define MIN_FLOAT (std::numeric_limits<float>::min())

class multinomial_regression_model_gradient_descent {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	std::size_t number_of_classes;
	float threshold; 							//the convergence for estimating the final weights
	float eta;									//used for gradient descent method
	MatrixXf experimental_results; 				//the experimental values of the features of the training data	
	MatrixXf experimental_results_trans;		//transpose of experimental_results
	MatrixXf weightsm; 							//weights of our learning network : F * K
	MatrixXf weightsm_trans;					//transpose of weights : K * F
	MatrixXf new_weightsm;						//updated weights after each step : F * K
	int* real_output;							//real output of each experimental results
	MatrixXf targets_multi_class;				//binary real output of each experimental results : N * K
	MatrixXf outputsm;							//outputs of the training data : N * K
	MatrixXf gradient;							//gradient of E : F * K
	MatrixXf sum_w_experimental_results;		//used for computing output : N 
	int* predicted_output_multi_class;			//predicted class of each experimental results
	

	void normalizing_weights_multi_class();
	void convert_target_to_binary(int* target_src, MatrixXf& targets_dst);	
	int eye_kj(std::size_t k, std::size_t j);
	void computing_all_output();
	void computing_all_gradient();
	void learning_weights_multi_classes();
	void new_values_for_weightsm();
	float computing_new_least_squared_err_multi_class();	
	void updating_values_of_weights_multi_class();
	void printing_weights_multi_class();
	void estimating_output_multiclass();
	
public:
	multinomial_regression_model_gradient_descent(std::size_t number_of_expr, std::size_t number_of_ftrs, std::size_t number_of_cls, 
											float th, float** expr_results, int* target_expr) {
		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs;
		number_of_classes = number_of_cls;
		threshold = th;
		eta = 0.001;	
	
		sum_w_experimental_results = MatrixXf::Random(number_of_experiments, 1);		
		weightsm = MatrixXf::Random(number_of_features, number_of_classes);
		weightsm_trans = MatrixXf::Random(number_of_classes, number_of_features);							
		new_weightsm = MatrixXf::Random(number_of_features, number_of_classes);
		gradient = MatrixXf::Random(number_of_features, number_of_classes);
		experimental_results = MatrixXf(number_of_experiments, number_of_features);
		experimental_results_trans = MatrixXf::Random(number_of_features, number_of_experiments);
		targets_multi_class = MatrixXf::Random(number_of_experiments, number_of_classes);
		outputsm = MatrixXf::Random(number_of_experiments, number_of_classes);
		predicted_output_multi_class = new int[number_of_experiments];
		real_output = new int[number_of_experiments];

		//initializing weights
		for(std::size_t f = 0; f < number_of_features; f++) {
			for(std::size_t k = 0; k < number_of_classes; k++) {
				weightsm(f, k) = 0.1;
			}
		}

		//initializing experimental_results
		// /experimental_results = expr_results;

		for(std::size_t i = 0; i < number_of_experiments; i++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
				experimental_results(i, f) = expr_results[i][f];
			}			
			//initializing real outputs
			real_output[i] = target_expr[i];
		}

		experimental_results_trans = experimental_results.transpose();
		
		//initializing targets_multi_class
		convert_target_to_binary(target_expr, targets_multi_class);
		outputsm = targets_multi_class;
	}

	void learning_multi_classes();
	MatrixXf& retrieving_weights_multi_classes();
	void printing_predicted_output_multi_class();
};

void multinomial_regression_model_gradient_descent::convert_target_to_binary(int* target_src, MatrixXf& targets_multi_class) {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t k = 0; k < number_of_classes; k++) {
			if(target_src[n] == k) {
				targets_multi_class(n, k) = 1.0;
			}
			else {
				targets_multi_class(n, k) = 0.0;
			}
		}
	}
}

//Ikj
int multinomial_regression_model_gradient_descent::eye_kj(std::size_t k, std::size_t j) {
	if(k == j) {
		return 1;
	}
	return 0;
}

//computing outputs
void multinomial_regression_model_gradient_descent::computing_all_output() {
	weightsm_trans = weightsm.transpose();
	//w^T * Q
	MatrixXf W_TQ_trans = MatrixXf::Random(number_of_classes, number_of_experiments);
	W_TQ_trans = weightsm_trans * experimental_results_trans;
	MatrixXf W_TQ = MatrixXf::Random(number_of_experiments, number_of_classes);
	W_TQ = W_TQ_trans.transpose();

	//sigma(exp(wQ))
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		sum_w_experimental_results(n, 0) = 0.0;
		for(std::size_t k = 0; k < number_of_classes; k++) {
			sum_w_experimental_results(n, 0) += exp(W_TQ(n, k));
		}
	}

	//ynk
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t k = 0; k < number_of_classes; k++) {		
			outputsm(n, k) = float(exp(W_TQ(n, k))/sum_w_experimental_results(n, 0)); 
		}
	}
}

//computing gradient 
void multinomial_regression_model_gradient_descent::computing_all_gradient(){
	//initializing
	gradient *= 0.0;
	for(std::size_t k = 0; k < number_of_classes; k++) {
		for(std::size_t n = 0; n < number_of_experiments; n++) {
			gradient.col(k) += (outputsm(n, k) - targets_multi_class(n, k)) * experimental_results_trans.col(n);
		}
	}
}

void multinomial_regression_model_gradient_descent::new_values_for_weightsm() {	
	new_weightsm = weightsm	- eta * gradient;
}

//computing leas squares err
float multinomial_regression_model_gradient_descent::computing_new_least_squared_err_multi_class() {	
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		if(predicted_output_multi_class[n] != real_output[n]) {
			num_err++;
		}		
	}
	float prec = float(num_err) / number_of_experiments;
	return prec;
}

//updating weights
void multinomial_regression_model_gradient_descent::updating_values_of_weights_multi_class() {	
	weightsm = new_weightsm;
}

void multinomial_regression_model_gradient_descent::printing_weights_multi_class() {
	for(std::size_t f = 0; f < number_of_features; f++) {
		for(std::size_t k = 0; k < number_of_classes; k++) {
			std::cout<<weightsm(f, k)<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n --------------------\n";
}

//estimating class of each experimental results based on the computed weights
void multinomial_regression_model_gradient_descent::estimating_output_multiclass() {	
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		float prob = MIN_FLOAT;
		for(std::size_t k = 0; k < number_of_classes; k++) {
			if(prob < outputsm(n, k)) {
				predicted_output_multi_class[n] = k;
				prob = outputsm(n, k);
			}
		}
	}
}

//updating weights till error meets the defined threshold
void multinomial_regression_model_gradient_descent::learning_weights_multi_classes() {
	float least_squared_err = MAX_FLOAT;
	std::size_t itr = 1;
	while(threshold < least_squared_err) {		
		computing_all_gradient();				
		new_values_for_weightsm();				
		updating_values_of_weights_multi_class();
		computing_all_output();
		estimating_output_multiclass();
		least_squared_err = computing_new_least_squared_err_multi_class();
		std::cout<<"("<<itr<<")"<<"Least_squared_err =\t" << least_squared_err<<std::endl;		
		printing_weights_multi_class();		
		itr++;
	}
}

void multinomial_regression_model_gradient_descent::normalizing_weights_multi_class() {
	float* averages = new float[number_of_features];
	float* averages_2 = new float[number_of_features];
	float* var = new float[number_of_features];
	//initializing
	for(std::size_t i = 0; i < number_of_features; i++) {
		averages[i] = 0;
		averages_2[i] = 0;
		var[i] = 0;
	}

	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features; j++) {
			averages[j] += experimental_results(i, j);
			averages_2[j] += (pow(experimental_results(i, j), 2.0));
		}
	}
	for(std::size_t i = 0; i < number_of_features; i++) {		
		averages[i] = float(averages[i]/number_of_experiments);
		averages_2[i] = float(averages_2[i]/number_of_experiments);
		var[i] = sqrt(averages_2[i] - pow(averages[i], 2.0));		
	}

	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
			experimental_results(n, f) = float((experimental_results(n, f) - averages[f])/var[f]);
		}
	}

	//releasing memory
	delete[] averages;
	delete[] averages_2;
	delete[] var;
}

void multinomial_regression_model_gradient_descent::learning_multi_classes() {
	normalizing_weights_multi_class();
	learning_weights_multi_classes();
}

MatrixXf& multinomial_regression_model_gradient_descent::retrieving_weights_multi_classes() {
	return weightsm;
}

void multinomial_regression_model_gradient_descent::printing_predicted_output_multi_class(){
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		std::cout<<"\n class["<<n<<"] =\t"<<predicted_output_multi_class[n]<<"\t"<<real_output[n];
		if(predicted_output_multi_class[n] != real_output[n]){
			num_err++;
		}
	}
	std::cout<<"\n number of error predicted is\t"<<num_err<<std::endl;
}
