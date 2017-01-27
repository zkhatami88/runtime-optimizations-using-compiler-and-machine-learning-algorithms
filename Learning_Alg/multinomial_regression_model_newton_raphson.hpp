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

class multinomial_regression_model_newton_raphson {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	std::size_t number_of_classes;
	float threshold; 							//the convergence for estimating the final weights
	MatrixXf experimental_results; 				//the experimental values of the features of the training data	
	MatrixXf experimental_results_trans;		//transpose of experimental_results
	MatrixXf weightsm; 							//weights of our learning network : F * K
	MatrixXf weightsm_trans;					//transpose of weights : K * F
	MatrixXf new_weightsm;						//updated weights after each step : F * K
	int* real_output;							//real output of each experimental results
	MatrixXf targets_multi_class;				//binary real output of each experimental results : N * K
	MatrixXf outputsm;							//outputs of the training data : N * K
	MatrixXf gradient;							//gradient of E : F * K
	MatrixXf hessian;							//hessian of E : (F * K) * (F * K)
	MatrixXf Ihessian;							//inverse of hessian of E : (F * K) * (F * K)
	MatrixXf delta_weight;						//used for updating weight : F * K
	MatrixXf sum_w_experimental_results;		//used for computing output : N 
	int* predicted_output_multi_class;			//predicted class of each experimental results
	

	void normalizing_weights_multi_class();
	void convert_target_to_binary(int* target_src, MatrixXf& targets_dst);	
	int eye_kj(std::size_t k, std::size_t j);
	void computing_all_output();
	void computing_all_gradient();
	void computing_heassian_kj(std::size_t k, std::size_t j);
	void computing_all_hessian();
	void learning_weights_multi_classes();
	void new_values_for_weightsm();
	float computing_new_least_squared_err_multi_class();	
	void updating_values_of_weights_multi_class();
	void printing_weights_multi_class();
	void estimating_output_multiclass();
	void printing_computed_values(std::size_t row, std::size_t col, MatrixXf& mat);
	
public:
	multinomial_regression_model_newton_raphson(std::size_t number_of_expr, std::size_t number_of_ftrs, std::size_t number_of_cls, 
											float th, float** expr_results, int* target_expr) {
		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs;
		number_of_classes = number_of_cls;
		threshold = th;	
	
		sum_w_experimental_results = MatrixXf::Random(number_of_experiments, 1);
		
		weightsm = MatrixXf::Random(number_of_features, number_of_classes);
		weightsm_trans = MatrixXf::Random(number_of_classes, number_of_features);							
		new_weightsm = MatrixXf::Random(number_of_features, number_of_classes);
		gradient = MatrixXf::Random(number_of_features, number_of_classes);
		delta_weight = 	MatrixXf::Random(number_of_features, number_of_classes);

		std::size_t num_row = number_of_features * number_of_classes;
		hessian = MatrixXf::Random(num_row, num_row);
		Ihessian = MatrixXf::Random(num_row, num_row);

		experimental_results = MatrixXf(number_of_experiments, number_of_features);
		experimental_results_trans = MatrixXf::Random(number_of_features, number_of_experiments);
		targets_multi_class = MatrixXf::Random(number_of_experiments, number_of_classes);
		outputsm = MatrixXf::Random(number_of_experiments, number_of_classes);
		predicted_output_multi_class = new int[number_of_experiments];
		real_output = new int[number_of_experiments];

		//initializing weights
		for(std::size_t f = 0; f < number_of_features; f++) {
			for(std::size_t k = 0; k < number_of_classes; k++) {
				weightsm(f, k) = 1;
			}
		}

		for(std::size_t n = 0; n < number_of_experiments; n++) {
			//initializing experimental_results
			for(std::size_t f = 0; f < number_of_features; f++) {
				experimental_results(n, f) = expr_results[n][f];
			}			
			//initializing real outputs
			real_output[n] = target_expr[n];
		}		
		
		//initializing targets_multi_class
		convert_target_to_binary(target_expr, targets_multi_class);
		outputsm = targets_multi_class;
	}

	void learning_multi_classes();
	void retrieving_weights_multi_classes_into_text_file();
	void printing_predicted_output_multi_class();
};

//it prints computed values : for testing
void multinomial_regression_model_newton_raphson::printing_computed_values(std::size_t row, std::size_t col, MatrixXf& mat) {
	if(row != 0 && col != 0) {
		for(std::size_t r = 0; r < row; r++) {
			for(std::size_t c = 0; c < col; c++) {
				printf("%f, ", mat(r, c));
			}
			std::cout<<std::endl;
		}
	}
	else if(row == 0 && col != 0){
		for(std::size_t c = 0; c < col; c++) {
			printf("%f, ", mat(0, c));
		}
	}
	else {
		for(std::size_t r = 0; r < row; r++) {
			printf("%f, ", mat(r, 0));
		}
	}
	std::cout<<std::endl;
}

void multinomial_regression_model_newton_raphson::convert_target_to_binary(int* target_src, MatrixXf& targets_multi_class) {
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
int multinomial_regression_model_newton_raphson::eye_kj(std::size_t k, std::size_t j) {
	if(k == j) {
		return 1;
	}
	return 0;
}

//computing outputs
void multinomial_regression_model_newton_raphson::computing_all_output() {
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
void multinomial_regression_model_newton_raphson::computing_all_gradient(){
	//initializing
	gradient *= 0.0;
	for(std::size_t k = 0; k < number_of_classes; k++) {
		for(std::size_t n = 0; n < number_of_experiments; n++) {
			gradient.col(k) += (outputsm(n, k) - targets_multi_class(n, k)) * experimental_results_trans.col(n);
		}
	}
}

void multinomial_regression_model_newton_raphson::computing_heassian_kj(std::size_t k, std::size_t j) {
	std::size_t row_offset = k * number_of_features;
	std::size_t col_offset = j * number_of_features;

	
	for(std::size_t r = 0; r < number_of_features; r++) {
		for(std::size_t c = 0; c < number_of_features; c++) {
			std::size_t row = row_offset + r;
			std::size_t col = col_offset + c;
			hessian(row, col) = 0.0;
			for(std::size_t n = 0; n < number_of_experiments; n++) {
				hessian(row, col) += outputsm(n, k) * (eye_kj(k, j) - outputsm(n, j)) * experimental_results(n, r) * experimental_results(n, c);
			}
			hessian(row, col) = (-1.0) * hessian(row, col);
		}
	}
}

void multinomial_regression_model_newton_raphson::computing_all_hessian() {
	for(std::size_t k = 0; k < number_of_classes; k++) {
		for(std::size_t j = 0; j < number_of_classes; j++) {
			computing_heassian_kj(k, j);
		}
	}
}

void multinomial_regression_model_newton_raphson::new_values_for_weightsm() {
	Ihessian = hessian.inverse();
	std::size_t size = number_of_classes * number_of_features;
	
	//computing delta_weight
	std::size_t r = 0;
	std::size_t counter_row = 0;
	std::size_t counter_col = 0;

	while(r < size) {
		delta_weight(counter_row, counter_col) = 0.0;		
		std::size_t h = 0;
		std::size_t which_class = 0;
		std::size_t counter = 0;
		while(h < size) {
			while(counter < number_of_features) {
					delta_weight(counter_row, counter_col) += hessian(r, h) * gradient(counter, which_class);
					counter++;
					h++;
				}
				counter = 0;
				which_class++;
		}
		counter_row++;
		r++;
		if(counter_row == number_of_features) {
			counter_row = 0;
			counter_col++;
		}
	}

	//updating weights
	new_weightsm = weightsm	- 1.2 * delta_weight;
}

//updating weights
void multinomial_regression_model_newton_raphson::updating_values_of_weights_multi_class() {	
	weightsm = new_weightsm;
}

void multinomial_regression_model_newton_raphson::printing_weights_multi_class() {
	printing_computed_values(number_of_features, number_of_classes, weightsm);
	std::cout<<"\n --------------------\n";
}

//estimating class of each experimental results based on the computed weights
void multinomial_regression_model_newton_raphson::estimating_output_multiclass() {	
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		float prob = outputsm(n, 0);
		for(std::size_t k = 1; k < number_of_classes; k++) {
			if(prob < outputsm(n, k)) {
				predicted_output_multi_class[n] = k;
				prob = outputsm(n, k);
			}
		}
	}
}

//computing leas squares err
float multinomial_regression_model_newton_raphson::computing_new_least_squared_err_multi_class() {	
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		if(predicted_output_multi_class[n] != real_output[n]) {
			num_err++;
		}		
	}
	float prec = float(num_err) / number_of_experiments;
	return prec;
}


//updating weights till error meets the defined threshold
void multinomial_regression_model_newton_raphson::learning_weights_multi_classes() {
	float least_squared_err = MAX_FLOAT;
	std::size_t itr = 1;
	while(threshold < least_squared_err) {
		computing_all_gradient();
		computing_all_hessian();
		new_values_for_weightsm();		
		updating_values_of_weights_multi_class();
		computing_all_output();
		estimating_output_multiclass();
		least_squared_err = computing_new_least_squared_err_multi_class();
		std::cout<<"("<<itr<<")"<<"Least_squared_err =\t" << least_squared_err<<std::endl;		
		printing_weights_multi_class();		
		itr++;
	}
	least_squared_err = computing_new_least_squared_err_multi_class();
	std::cout<<"("<<itr<<") : "<<"Least_squared_err =\t" << least_squared_err<<std::endl;
}

void multinomial_regression_model_newton_raphson::normalizing_weights_multi_class() {
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

void multinomial_regression_model_newton_raphson::learning_multi_classes() {
	normalizing_weights_multi_class();
	experimental_results_trans = experimental_results.transpose();
	learning_weights_multi_classes();
}

void multinomial_regression_model_newton_raphson::retrieving_weights_multi_classes_into_text_file() {
	std::ofstream outputFile("weights.txt");
	for(std::size_t f = 0; f < number_of_features; f++) {
		for(std::size_t k = 0; k < number_of_classes - 1; k++) {
			outputFile<<weightsm(f, k)<<" ";
		}
		outputFile<<weightsm(f, number_of_classes - 1)<<" ";
		if(f != number_of_features - 1) {
			outputFile<<std::endl;
		}
	}
}

void multinomial_regression_model_newton_raphson::printing_predicted_output_multi_class(){
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		if(predicted_output_multi_class[n] != real_output[n]){
			num_err++;
			std::cout<<"\n["<<n<<"] =\t"<<predicted_output_multi_class[n]<<"\t"<<real_output[n];
		}
	}
	std::cout<<"\n number of error predicted is\t"<<num_err<<" out of "<<number_of_experiments<<std::endl;
}
