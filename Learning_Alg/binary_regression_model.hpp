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

class learning_binary_regression_model {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	float threshold; 							//the convergence for estimating the final weights
	MatrixXf experimental_results; 				//the experimental values of the features of the training data	
	MatrixXf experimental_results_trans;		//transpose of experimental_results
	MatrixXf weightsb; 							//weights of our learning network
	MatrixXf new_weightsb;						//updated weights after each step	
	MatrixXf targets_two_class;					//outputs of the training data
	MatrixXf diag_weightsb; 					//used for updating weights
	MatrixXf M; 								//used for updating weights
	int* predicted_output_two_class;			//predicted class of each experimental results
		
	
	void normalizing_samples_two_class();
	void updating_values_of_M_and_diag_weights();
	void updating_values_of_weights_two_class();
	void new_values_for_weightsb();
	float computing_new_least_squared_err_two_class();	
	void learning_weights_two_classes();
	void printing_weights_two_class();
	void estimating_output_two_class();
	void printing_computed_values(std::size_t row, std::size_t col, MatrixXf& mat);

public:
	learning_binary_regression_model(std::size_t number_of_expr, std::size_t number_of_ftrs, 
									float th, float** expr_results, int* target_expr) {
		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs;
		threshold = th;

		weightsb = Eigen::MatrixXf::Random((number_of_features +  1), 1);
		new_weightsb = Eigen::MatrixXf::Random((number_of_features +  1), 1);
		targets_two_class = Eigen::MatrixXf::Random(number_of_experiments, 1);
		diag_weightsb = Eigen::MatrixXf::Random(number_of_experiments, number_of_experiments);
		M = Eigen::MatrixXf::Random(number_of_experiments, 1);
		experimental_results = Eigen::MatrixXf::Random(number_of_experiments, (number_of_features +  1));
		experimental_results_trans = Eigen::MatrixXf::Random((number_of_features + 1), number_of_experiments);	
		predicted_output_two_class = new int[number_of_experiments];

		//initializing weights
		for(std::size_t i = 0; i < number_of_features +  1; i++) {
			weightsb(i, 0) = 0.1;
		}

		for(std::size_t n = 0; n < number_of_experiments; n++) {
			targets_two_class(n, 0) = (float)target_expr[n];
			M(n, 0) = 0.0;

			//initializing experimental_results
			experimental_results(n, 0) = 1.0; // 1 + wf + ....
			for(std::size_t f = 1; f < number_of_features + 1; f++) {
				experimental_results(n, f) = expr_results[n][f - 1];
			}

			//initializing diag_weightsb
			for(std::size_t j = 0; j < number_of_experiments; j++) {
				diag_weightsb(n, j) = 0.0;
			}
		}		
	}

	//two-class logistic regression
	void learning_two_classes();
	MatrixXf& retrieving_weights_two_classes();
	void finalizing_two_classes();
	void printing_predicted_output_two_class();
};

//it prints computed values : for testing
void learning_binary_regression_model::printing_computed_values(std::size_t row, std::size_t col, MatrixXf& mat) {
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

void learning_binary_regression_model::updating_values_of_M_and_diag_weights() {
	MatrixXf WX = MatrixXf::Random(1, number_of_experiments);
	MatrixXf weightsb_transpose = MatrixXf::Random(1, (number_of_features + 1));
	WX = weightsb_transpose * experimental_results_trans;

	//std::cout<<"\n ---experimental_results:\n";
	//printing_computed_values(number_of_experiments, number_of_features + 1, experimental_results);

	//std::cout<<"\n ---WX : \n";
	//printing_computed_values(0, number_of_experiments, WX);

	for(std::size_t n = 0; n < number_of_experiments; n++) {
		M(n, 0) = float(1.0 / (1.0 + exp((-1.0) *WX(0, n))));
		diag_weightsb(n, n) = M(n, 0) * (1.0 - M(n, 0));
	}

	//std::cout<<"\n ----M : \n";
	//printing_computed_values(number_of_experiments, 0, M);

	//std::cout<<"\n ----diag_weightsb : \n";
	//printing_computed_values(number_of_experiments, number_of_experiments, diag_weightsb);
}

void learning_binary_regression_model::new_values_for_weightsb() {
	//X^T * S
	MatrixXf X_TS = Eigen::MatrixXf::Random((number_of_features + 1), number_of_experiments);
	X_TS = experimental_results_trans * diag_weightsb;

	//std::cout<<"\n ----X_TS : \n";
	//printing_computed_values(number_of_features + 1, number_of_experiments, X_TS);

	//X^T * S * X
	MatrixXf XSX = Eigen::MatrixXf::Random((number_of_features + 1), (number_of_features + 1));
	XSX = X_TS * experimental_results;

	//std::cout<<"\n ----XSX : \n";
	//printing_computed_values(number_of_features + 1, number_of_features + 1, XSX);

	//(X^T * S * X) ^ -1
	MatrixXf XSX_inv = Eigen::MatrixXf::Random((number_of_features + 1), (number_of_features + 1));
	XSX_inv = XSX.inverse();

	//std::cout<<"\n ----XSX_inv : \n";
	//printing_computed_values(number_of_features + 1, number_of_features + 1, XSX_inv);

	//(X^T * S * X) ^ -1  * X^T
	MatrixXf XSX_inv_X_T = Eigen::MatrixXf::Random((number_of_features + 1), number_of_experiments);
	XSX_inv_X_T = XSX_inv * experimental_results_trans;

	//std::cout<<"\n ----XSX_inv_X_T : \n";
	//printing_computed_values(number_of_features + 1, number_of_experiments, XSX_inv_X_T);

	//S * X
	MatrixXf SX = Eigen::MatrixXf::Random(number_of_experiments, (number_of_features + 1));
	SX = diag_weightsb * experimental_results;

	//std::cout<<"\n ----SX : \n";
	//printing_computed_values(number_of_experiments, number_of_features + 1, SX);

	//S * X * W
	MatrixXf SXW = Eigen::MatrixXf::Random(number_of_experiments, 1);
	SXW = SX * weightsb;
	
	//std::cout<<"\n ----SXW : \n";
	//printing_computed_values(number_of_experiments, 0, SXW);

	//S * X * W + y - M
	MatrixXf SXW_y_M = Eigen::MatrixXf::Random(number_of_experiments, 1);
	SXW_y_M = SXW + targets_two_class - M;

	//std::cout<<"\n ----SXW_y_M : \n";
	//printing_computed_values(number_of_experiments, 0, SXW_y_M);

	//W
	new_weightsb = XSX_inv_X_T * SXW_y_M;

	//std::cout<<"\n ----new_weightsb : \n";
	//printing_computed_values(number_of_features + 1, 0, new_weightsb);
}

float learning_binary_regression_model::computing_new_least_squared_err_two_class() {
	float num_err = 0.0;
	for(int n = 0; n < number_of_experiments; n++) {
		if(predicted_output_two_class[n] != targets_two_class(n, 0)) {
			num_err++;
		}
	}
	float prec = float(num_err / number_of_experiments);
	return prec;
}

void learning_binary_regression_model::updating_values_of_weights_two_class() {
	for(std::size_t f = 0; f < number_of_features + 1; f++) {
		weightsb(f, 0) = new_weightsb(f, 0);
	}
}

void learning_binary_regression_model::printing_weights_two_class() {
	printing_computed_values(number_of_features + 1, 0, weightsb);
}

void learning_binary_regression_model::estimating_output_two_class(){
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		float temp = 0.0;
		for(std::size_t f = 0; f < number_of_features + 1; f++) {
			temp += weightsb(f, 0) * experimental_results(n, f);
		}
		if(temp >= 0) {
			predicted_output_two_class[n] = 1;
		}
		else {
			predicted_output_two_class[n] = 0;
		}
	}
}

//this func applys experimental values of the training data on our learning network
void learning_binary_regression_model::learning_weights_two_classes() {
	float least_squared_err = MAX_FLOAT;	
	std::size_t itr = 1;
	while(threshold < least_squared_err) {	
		updating_values_of_M_and_diag_weights();
		new_values_for_weightsb();
		updating_values_of_weights_two_class();	
		estimating_output_two_class();
		least_squared_err = computing_new_least_squared_err_two_class();
		std::cout<<"\n("<<itr<<"),"<<"Least_squared_err =\t" << least_squared_err<<std::endl;			
		printing_weights_two_class();		
		itr++;		
	}	
}

void learning_binary_regression_model::normalizing_samples_two_class() {
	//normalizing all features, except first one, which is 1 for all
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
		for(std::size_t j = 1; j < number_of_features + 1; j++) {
			averages[j - 1] += experimental_results(i, j);
			averages_2[j - 1] += (pow(experimental_results(i, j), 2.0));
		}
	}
	for(std::size_t i = 0; i < number_of_features; i++) {		
		averages[i] = float(averages[i]/number_of_experiments);
		averages_2[i] = float(averages_2[i]/number_of_experiments);
		var[i] = sqrt(averages_2[i] - pow(averages[i], 2.0));		
	}

	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t f = 1; f < number_of_features + 1; f++) {
			if(var[f - 1] != 0) {
				experimental_results(n, f) = float((experimental_results(n, f) - averages[f - 1])/var[f - 1]);
			}
		}
	}

	//releasing memory	
	delete[] averages;
	delete[] averages_2;
	delete[] var;
}

void learning_binary_regression_model::learning_two_classes() {
	normalizing_samples_two_class();
	experimental_results_trans = experimental_results.transpose();
	learning_weights_two_classes();
}

MatrixXf& learning_binary_regression_model::retrieving_weights_two_classes() {
	printing_weights_two_class();
	return weightsb;
}

void learning_binary_regression_model::printing_predicted_output_two_class() {
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		if(predicted_output_two_class[n] != targets_two_class(n, 0)) {
			num_err++;
			std::cout<<"["<<n<<"]\t"<<predicted_output_two_class[n]<<","<<targets_two_class(n, 0)<<std::endl;
		}		
	}
	std::cout<<"\n number of error predicted is\t"<<num_err<<" out of "<<number_of_experiments<<std::endl;
}
