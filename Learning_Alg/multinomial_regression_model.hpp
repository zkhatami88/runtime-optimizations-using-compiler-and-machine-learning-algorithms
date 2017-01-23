#include <limits>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h> 
/*
computing inverse of matrix
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/LU>

using namespace Eigen;


#define MAX_double (std::numeric_limits<double>::max())
#define MIN_double (std::numeric_limits<double>::min())

class learning_network {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	std::size_t number_of_classes;
	double threshold; 							//the convergence for estimating the final weights
	double** experimental_results; 				//the experimental values of the features of the training data
	

	//two-class logistic regression
	double* weightsb; 							//weights of our learning network
	double* new_weightsb;						//updated weights after each step	
	double* targets_two_class;				//outputs of the training data
	double** diag_weightsb; 					//used for updating weights
	double* M; 									//used for updating weights
	int* predicted_output_two_class;					//predicted class of each experimental results
	
	//multi-class logistic regression
	double** weightsm; 							//weights of our learning network : (F * K) * 1
	double** new_weightsm;						//updated weights after each step : (F * K) * 1
	int** targets_multi_class;				//real output of each experimental results : N * K
	double** outputsm;							//outputs of the training data : N * K
	double** gradient;							//gradient of E : (F * K) * 1
	double** hessian;							//hessian of E : (F * K) * (F * K)
	double** Ihessian;							//inverse of hessian of E : (F * K) * (F * K)
	double** delta_weight;						//used for updating weight : (F * K) * 1
	double* sum_w_experimental_results;			//used for computing output : N 
	int* predicted_output_multi_class;					//predicted class of each experimental results
	int* real_output;						//real output of each experimental results
	
	void inverse_mat(double** mat, double** inverse, std::size_t size);
	void get_cofactor(double **src, double **dest, std::size_t p, std::size_t q, std::size_t size);
	double calc_determinant( double **mat, std::size_t size);
	void adjoint(double** mat, double** adj, std::size_t size);
	
	/*
	two-class logistic regression (binary logistic regression) :
	*/
	void normalizing_weights_two_class();
	double get_value_of_M(std::size_t i);
	void updating_values_of_M_and_diag_weights();		
	void inverse_martix_XSX(double** mat, double** inverse_mat, std::size_t size);
	void updating_values_of_weights_two_class();
	void new_values_for_weightsb();
	double computing_new_least_squared_err_two_class();	
	void learning_weights_two_classes();
	void printing_weights_two_class();
	void estimating_output_two_class();

	/*
	multi-class logistic regression 
	*/
	void normalizing_weights_multi_class();
	void convert_target_to_binary(int* target_src, int** targets_dst);
	void sum_w_experimental_results_n(std::size_t n);
	double output_nk(std::size_t n, std::size_t k);
	void computing_all_output_n(std::size_t n);
	void computing_all_output();
	std::size_t eye_kj(std::size_t k, std::size_t j);
	void gradient_E_k(std::size_t k);
	void computing_all_gradient();
	double hessian_E_kj(std::size_t k, std::size_t j);
	void computing_all_hessian();
	void eigen_mat_initilializer(MatrixXf&& eig_mat, std::size_t size);
	void set_Ihessian_value(MatrixXf&& inv_eig_mat, std::size_t size);
	void computing_inverse_hessian();
	void learning_weights_multi_classes();
	void computing_delta_weight();
	void new_values_for_weightsm();
	double computing_new_least_squared_err_multi_class();	
	void updating_values_of_weights_multi_class();
	void printing_weights_multi_class();
	void estimating_output_multiclass();
	
public:
	learning_network(std::size_t number_of_expr, std::size_t number_of_ftrs, 
						std::size_t number_of_cls, double th) {
		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs;
		number_of_classes = number_of_cls;
		threshold = th;		
	}	

	//initilializer for two_class
	void initilializer_two_class(double** expr_results, int* target_expr){
		weightsb = new double[number_of_features +  1];
		new_weightsb = new double[number_of_features +  1];
		targets_two_class = new double[number_of_experiments];
		diag_weightsb = new double*[number_of_experiments];
		M = new double[number_of_experiments];
		experimental_results = new double*[number_of_experiments];
		predicted_output_two_class = new int[number_of_experiments];

		//initializing weights
		for(std::size_t i = 0; i < number_of_features +  1; i++) {
			weightsb[i] = 0.1;
		}				

		for(std::size_t i = 0; i < number_of_experiments; i++) {
			experimental_results[i] = new double[number_of_features + 1];
			targets_two_class[i] = (double)target_expr[i];
			diag_weightsb[i] = new double[number_of_experiments];
			M[i] = 0.0;

			//initializing experimental_results
			experimental_results[i][0] = 1.0; // 1 + wf + ....
			for(std::size_t j = 1; j < number_of_features + 1; j++) {
				experimental_results[i][j] = expr_results[i][j - 1];
			}

			//initializing diag_weightsb
			for(std::size_t j = 0; j < number_of_experiments; j++) {
				diag_weightsb[i][j] = 0.0;
			}
		}

		for(std::size_t i = 0; i < number_of_experiments; i++) {
			for(std::size_t j = 0; j < number_of_features + 1; j++) {
				std::cout<<experimental_results[i][j]<<"\t";
			}
			std::cout<<targets_two_class[i];
			std::cout<<std::endl;
		}
	}	

	//initilializer for multi_class
	void initilializer_multi_class(double** expr_results, int* target_expr) {
		sum_w_experimental_results = new double[number_of_experiments];
		std::size_t num_row = number_of_features * number_of_classes;
		weightsm = new double*[num_row]; 							
		new_weightsm = new double*[num_row];
		gradient = new double*[num_row];
		delta_weight = 	new double*[num_row];
		hessian = new double*[num_row];
		Ihessian = new double*[num_row];
		experimental_results = new double*[number_of_experiments];
		targets_multi_class = new int*[number_of_experiments];
		outputsm = new double*[number_of_experiments];
		predicted_output_multi_class = new int[number_of_experiments];
		real_output = new int[number_of_experiments];

		for(std::size_t i = 0; i < num_row; i++) {
			weightsm[i] =  new double[1];
			new_weightsm[i] =  new double[1];
			gradient[i] =  new double[1];
			delta_weight[i] =  new double[1];
			hessian[i] =  new double[num_row];
			Ihessian[i] =  new double[num_row];

			//initializing weights
			srand (i);
			weightsm[i][0] =  (double)(rand() % 10)/10 + 0.1;
		}
		
		for(std::size_t i = 0; i < number_of_experiments; i++) {
			experimental_results[i] = new double[number_of_features];
			targets_multi_class[i] = new int[number_of_classes];
			outputsm[i] = new double[number_of_classes];

			//initializing experimental_results
			for(std::size_t f = 0; f < number_of_features; f++) {
				experimental_results[i][f] = expr_results[i][f];
			}

			//initializing real outputs
			real_output[i] = target_expr[i];
		}

		//initializing targets_multi_class
		convert_target_to_binary(target_expr, targets_multi_class);
	}

	//two-class logistic regression
	void learning_two_classes();
	double* retrieving_weights_two_classes();
	void finalizing_two_classes();
	void printing_predicted_output_two_class();

	//multi-class logistic regression
	void learning_multi_classes();
	double** retrieving_weights_multi_classes();
	void finalizing_multi_classes();
	void printing_predicted_output_multi_class();
};


/*
two-class logistic regression (binary logistic regression) :
*/

double learning_network::get_value_of_M(std::size_t i) {
	double temp = 0.0;
	for(std::size_t j = 0; j < number_of_features + 1; j++) {
		temp += weightsb[j] * experimental_results[i][j];
	}
	temp = temp * (-1);
	double result = double(1 / (1 + exp(temp)));
	return result;
}

void learning_network::updating_values_of_M_and_diag_weights() {
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		M[i] = get_value_of_M(i);
		diag_weightsb[i][i] = M[i] * (1 - M[i]);
	}
}


void learning_network::inverse_martix_XSX(double** mat, double** inverse_mat, std::size_t size) {
	MatrixXf eig_mat = MatrixXf::Random(size, size);
	MatrixXf inv_eig_mat = MatrixXf::Random(size, size);

	std::cout<<"\n SIZE = "<<size<<std::endl;
	std::cout<<"\n------------------------------\n";
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j  = 0; j <size; j++) {
			std::cout<<mat[i][j]<<",";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n------------------------------\n";

	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			eig_mat(i, j) = mat[i][j];
		}
	}

	inv_eig_mat = eig_mat.inverse();

	std::cout<<"\n----------inverse----------\n";
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j  = 0; j <size; j++) {
			std::cout<<inv_eig_mat(i, j)<<",";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n------------------------------\n";

	for(std::size_t i = 0 ; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			inverse_mat[i][j] = inv_eig_mat(i ,j); 
		}
	}

    //inverse_mat(mat, inverse, size);
}

void learning_network::new_values_for_weightsb() {
	//allocating memory
	double** temp4; //X^T * S
	double** temp1; //XSX : f * f
	double** temp2; //SXW + y - M : expr * 1
	double** temp3; //XSX * X^T : f * expr	
	double** inverse_temp1;
	temp1 = new double*[number_of_features + 1];
	temp2 = new double*[number_of_experiments];
	temp3 = new double*[number_of_features + 1];
	temp4 = new double*[number_of_features + 1]; 
	inverse_temp1 = new double*[number_of_features + 1];    

	//initializing
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
    	inverse_temp1[i] = new double[number_of_features + 1];
    }

	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		temp1[i] = new double[number_of_features + 1];
		temp3[i] = new double[number_of_experiments];
		temp4[i] = new double[number_of_experiments];
	}	
	
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		temp2[i] = new double[1];
		temp2[i][0] = 0.0;
	}

	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		for(std::size_t j = 0; j < number_of_experiments; j++) {			
			temp4[i][j] = experimental_results[j][i] * diag_weightsb[j][j];
		}

		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			for(std::size_t k = 0; k < number_of_experiments; k++) {
				//initializing value
				if(k == 0) {
					temp1[i][j] = 0;
				}
				temp1[i][j] += temp4[i][k] * experimental_results[k][j];
			}
		}
	}

	inverse_martix_XSX(temp1, inverse_temp1, number_of_features + 1);

	//computing temp2
	for(std::size_t i = 0; i < number_of_experiments; i++) { 
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			temp2[i][0] += diag_weightsb[i][i] * experimental_results[i][j] * weightsb[j]; 
		}
		temp2[i][0] += targets_two_class[i] - M[i];
	}

	//updating weightsb
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		new_weightsb[i] = 0.0;
		for(std::size_t j = 0; j < number_of_experiments; j++) {
			for(std::size_t k = 0; k < number_of_features + 1; k++) {
				if(k == 0) {
					temp3[i][j] = 0.0;
				}
				temp3[i][j] += inverse_temp1[i][k] * experimental_results[j][k];
			}
			new_weightsb[i] += temp3[i][j] * temp2[j][0];
		}
	}

	//releasing memory
	/*
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		delete[] temp1[i];
		delete[] inverse_temp1[i];
		delete[] temp3[i];
		delete[] temp4[i];
		temp1[i] = nullptr;
		inverse_temp1[i] = nullptr;
		temp3[i] = nullptr;
		temp4[i] = nullptr;
	}
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		delete[] temp2[i];
		temp2[i] = nullptr;
	}
	delete[] temp1;
	delete[] inverse_temp1;
	delete[] temp2;
	delete[] temp3;
	delete[] temp4;*/
}

double learning_network::computing_new_least_squared_err_two_class() {
	double temp_err = 0.0;
	/*
	double old_wx, new_wz;
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		old_wx = 0;
		new_wz = 0;
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			old_wx += weightsb[j] * experimental_results[i][j];
			new_wz += new_weightsb[j] * experimental_results[i][j];
		}
		temp_err += pow( ((1 / (1 + (1 / exp(old_wx))))) - (1 / (1 + (1 / exp(new_wz)))) , 2.0);
	}*/
	for(std::size_t f = 0; f < number_of_features + 1; f++) {
		temp_err += pow((weightsb[f] - new_weightsb[f]), 2.0);
	}
	return sqrt(temp_err);
}

void learning_network::updating_values_of_weights_two_class() {
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		weightsb[i] = new_weightsb[i];
	}
}

void learning_network::printing_weights_two_class() {
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		std::cout<<"weights["<<i<<"] = "<<weightsb[i]<<"\t";
	}
	std::cout<<"\n----------"<<std::endl;
}

//this func applys experimental values of the training data on our learning network
void learning_network::learning_weights_two_classes() {
	double least_squared_err = MAX_double;	
	std::size_t itr = 0;
	//while(threshold < least_squared_err) {		
		updating_values_of_M_and_diag_weights();
		new_values_for_weightsb();
		least_squared_err = computing_new_least_squared_err_two_class();
		std::cout<<"("<<itr<<")\t"<<"Least_squared_err =\t" << least_squared_err<<std::endl;
		updating_values_of_weights_two_class();		
		printing_weights_two_class();
		itr++;
	//}
}

void learning_network::normalizing_weights_two_class() {
	//normalizing all features, except first one, which is 1 for all
	double* averages = new double[number_of_features];
	double* averages_2 = new double[number_of_features];
	double* var = new double[number_of_features];
	//initializing
	for(std::size_t i = 0; i < number_of_features; i++) {
		averages[i] = 0;
		averages_2[i] = 0;
		var[i] = 0;
	}

	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 1; j < number_of_features + 1; j++) {
			averages[j - 1] += experimental_results[i][j];
			averages_2[j - 1] += (pow(experimental_results[i][j], 2.0));
		}
	}
	for(std::size_t i = 0; i < number_of_features; i++) {		
		averages[i] = double(averages[i]/number_of_experiments);
		averages_2[i] = double(averages_2[i]/number_of_experiments);
		var[i] = sqrt(averages_2[i] - pow(averages[i], 2.0));		
	}

	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t f = 1; f < number_of_features + 1; f++) {
			experimental_results[n][f] = (experimental_results[n][f] - averages[f])/var[f];
		}
	}

	//releasing memory
	/*
	delete[] averages;
	delete[] averages_2;
	delete[] var;*/
}

void learning_network::learning_two_classes() {
	//normalizing_weights_two_class();
	learning_weights_two_classes();
}

double* learning_network::retrieving_weights_two_classes() {
	return weightsb;
}

void learning_network::estimating_output_two_class(){
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		double temp = 0.0;
		for(std::size_t f = 0; f < number_of_features + 1; f++) {
			temp += weightsb[f] * experimental_results[n][f];
		}
		if(temp > 0) {
			predicted_output_two_class[n] = 0;
		}
		else {
			predicted_output_two_class[n] = 1;
		}
	}
}

void learning_network::printing_predicted_output_two_class() {
	estimating_output_two_class();
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		std::cout<<predicted_output_two_class[n]<<"\t"<<targets_two_class[n]<<std::endl;
	}
}

void learning_network::finalizing_two_classes() {
	//releasing memory
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		delete[] experimental_results[n];
		delete[] diag_weightsb[n];
		experimental_results[n] = nullptr;
		diag_weightsb[n] = nullptr;
	}		
	delete[] experimental_results;
	delete[] diag_weightsb;
	delete[] weightsb;
	delete[] new_weightsb;
	delete[] targets_two_class;		
	delete[] M;
	delete[] predicted_output_two_class;
	weightsb = nullptr;
	new_weightsb = nullptr;
	targets_two_class = nullptr;
	M = nullptr;
	predicted_output_two_class = nullptr;
}

/*
multi-class logistic regression 
*/

void learning_network::convert_target_to_binary(int* target_src, int** targets_dst) {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t k = 0; k < number_of_classes; k++) {
			if(target_src[n] == k) {
				targets_dst[n][k] = 1;
			}
			else {
				targets_dst[n][k] = 0;
			}
		}
	}
}

//Ikj
std::size_t learning_network::eye_kj(std::size_t k, std::size_t j) {
	if(k == j) {
		return 1;
	}
	return 0;
}

//computing ouput
void learning_network::sum_w_experimental_results_n(std::size_t n) {
	sum_w_experimental_results[n] = 0.0;	
	for(std::size_t k = 0; k < number_of_classes; k++) {
		std::size_t offset = k * number_of_features;
		double temp_sum = 0.0;
		for(std::size_t f = 0; f < number_of_features; f++) {
			temp_sum += weightsm[offset + f][0] * experimental_results[n][f];
		}
		sum_w_experimental_results[n] += exp(temp_sum);
	}
}

double learning_network::output_nk(std::size_t n, std::size_t k) {
	double temp = 0.0;
	std::size_t offset = k * number_of_features;
	for(std::size_t f = 0; f < number_of_features; f++) {
		temp += weightsm[offset + f][0] * experimental_results[n][f];
	}

	double out = (double)((double)exp(temp)/(double)sum_w_experimental_results[n]);
	if(n == 0){
		std::cout<<"\n sum_w_experimental_results["<<n<<"] ="<<sum_w_experimental_results[n]<<", outputsm =\t"<<out<<std::endl;
	}
	return out;
}

void learning_network::computing_all_output_n(std::size_t n) {	
	sum_w_experimental_results_n(n);
	for(std::size_t k = 0; k < number_of_classes; k++) {
		outputsm[n][k] = output_nk(n, k);
	}
}

void learning_network::computing_all_output() {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		computing_all_output_n(n);		
	}
}

//computing gradient for each W
void learning_network::gradient_E_k(std::size_t k) {
	std::size_t offset = k * number_of_features;
	for(std::size_t f = 0; f < number_of_features; f++) {
		gradient[offset + f][0] = 0.0;
		for(std::size_t n = 0; n < number_of_experiments; n++) {
			gradient[offset + f][0] += (output_nk(n, k) - (double)targets_multi_class[n][k]) * experimental_results[n][f];			
		}
	}
}

void learning_network::computing_all_gradient(){
	for(std::size_t k = 0; k < number_of_classes; k++) {
		gradient_E_k(k);
	}
}

//computing hessian 
double learning_network::hessian_E_kj(std::size_t k, std::size_t j) {
	std::size_t row_offset = k * number_of_features;
	std::size_t col_offset = j * number_of_features;

	for(std::size_t i = 0; i < number_of_features; i++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
			hessian[row_offset + i][col_offset + f] = 0.0;
			for(std::size_t n = 0; n < number_of_experiments; n++) {
				hessian[row_offset + i][col_offset + f] += outputsm[n][k] * (eye_kj(k, j) - outputsm[n][j]) * experimental_results[n][i] * experimental_results[n][f];
			}
			hessian[row_offset + i][col_offset + f] = (-1) * hessian[row_offset + i][col_offset + f];
		}
	}
}

void learning_network::computing_all_hessian() {
	for(std::size_t k  = 0; k < number_of_classes; k++) {
		for(std::size_t j = 0; j < number_of_classes; j++) {
			hessian_E_kj(k, j);
		}
	}
	std::cout<<"\n ================================= \n";
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			std::cout<<hessian[i][j]<<",";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n ================================= \n";
}

//initializing eig_mat with hessian values
void learning_network::eigen_mat_initilializer(MatrixXf&& eig_mat, std::size_t size) {
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			eig_mat(i, j) = hessian[i][j];
		}
	}
}

//setting Ihessian values with inverse of eig_mat, which is hessian
void learning_network::set_Ihessian_value(MatrixXf&& inv_eig_mat, std::size_t size) {
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			Ihessian[i][j] = inv_eig_mat(i, j);
		}
	}
}

void learning_network::computing_inverse_hessian() {
	std::size_t size = number_of_classes * number_of_features;
	MatrixXf eig_mat = MatrixXf::Random(size, size);
	MatrixXf inv_eig_mat = MatrixXf::Random(size, size);
	eigen_mat_initilializer(std::move(eig_mat), size);
	inv_eig_mat = eig_mat.inverse();
	set_Ihessian_value(std::move(inv_eig_mat), size);
}

//Ihessian * gradient
void learning_network::computing_delta_weight() {
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t i = 0; i < size; i++) {
		delta_weight[i][0] = 0;
		for(std::size_t j = 0; j < size; j++) {
			delta_weight[i][0] += Ihessian[i][j] * gradient[j][0];
		}
	}
}

void learning_network::new_values_for_weightsm() {
	computing_inverse_hessian();
	std::size_t num_row = number_of_classes * number_of_features;
	std::cout<<"\n ----------------Inverse---------------\n";
	for(std::size_t i = 0; i < num_row; i++) {
		for(std::size_t j = 0; j < num_row; j++){
			std::cout<<Ihessian[i][j]<<",";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n ----------------------------- \n";

	computing_delta_weight();
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t i = 0; i < size; i++) {
		new_weightsm[i][0] = weightsm[i][0] - delta_weight[i][0];
	}
}

//computing leas squares err
double learning_network::computing_new_least_squared_err_multi_class() {
	double err = 0.0;
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t i = 0; i < size; i++) {
		err += pow((new_weightsm[i][0] - weightsm[i][0]), 2.0);
	}
	return sqrt(err);
}

//updating weights
void learning_network::updating_values_of_weights_multi_class() {
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t i = 0; i < size; i++) {
		weightsm[i][0] = new_weightsm[i][0];
	}
}

void learning_network::printing_weights_multi_class() {
	std::size_t size = number_of_classes * number_of_features;
	for(std::size_t k = 0; k < number_of_classes; k++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
			std::size_t offset = k * number_of_features + f;
			//std::cout<<"weights["<<k<<"]["<<f<<"] = "<<weightsm[offset][0]<<"\t";
			//std::cout<<"gradient["<<k<<"]["<<f<<"] = "<<gradient[offset][0]<<"\t";
			//std::cout<<"delta_weight["<<k<<"]["<<f<<"] = "<<delta_weight[offset][0]<<"\t";
			//std::cout<<std::endl;
		}
				
	}
}

//updating weights till error meets the defined threshold
void learning_network::learning_weights_multi_classes() {
	double least_squared_err = MAX_double;
	std::size_t itr = 0;
	while(threshold < least_squared_err) {
		computing_all_output();		
		computing_all_gradient();
		computing_all_hessian();				
		new_values_for_weightsm();				
		least_squared_err = computing_new_least_squared_err_multi_class();
		std::cout<<"("<<itr<<")"<<"Least_squared_err =\t" << least_squared_err<<std::endl;		
		updating_values_of_weights_multi_class();
		printing_weights_multi_class();		
		itr++;
	}
}

void learning_network::normalizing_weights_multi_class() {
	double* averages = new double[number_of_features];
	double* averages_2 = new double[number_of_features];
	double* var = new double[number_of_features];
	//initializing
	for(std::size_t i = 0; i < number_of_features; i++) {
		averages[i] = 0;
		averages_2[i] = 0;
		var[i] = 0;
	}

	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features; j++) {
			averages[j] += experimental_results[i][j];
			averages_2[j] += (pow(experimental_results[i][j], 2.0));
		}
	}
	for(std::size_t i = 0; i < number_of_features; i++) {		
		averages[i] = double(averages[i]/number_of_experiments);
		averages_2[i] = double(averages_2[i]/number_of_experiments);
		var[i] = sqrt(averages_2[i] - pow(averages[i], 2.0));		
	}

	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
			experimental_results[n][f] = (experimental_results[n][f] - averages[f])/var[f];
		}
	}

	//releasing memory
	delete[] averages;
	delete[] averages_2;
	delete[] var;
}

void learning_network::learning_multi_classes() {
	//normalizing_weights_multi_class();
	learning_weights_multi_classes();
}

double** learning_network::retrieving_weights_multi_classes() {
	return weightsm;
}

//estimating class of each experimental results based on the computed weights
void learning_network::estimating_output_multiclass() {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		double prob = MIN_double;
		for(std::size_t k = 0; k < number_of_classes; k++) {
			double out = output_nk(n, k);
			if(prob < outputsm[n][k]) {
				predicted_output_multi_class[n] = k;
				prob = outputsm[n][k];
			}
		}
	}
}

void learning_network::printing_predicted_output_multi_class(){
	estimating_output_multiclass();
	std::size_t num_err = 0;
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		std::cout<<"\n class["<<n<<"] =\t"<<predicted_output_multi_class[n]<<"\t"<<real_output[n];
		if(predicted_output_multi_class[n] != real_output[n]){
			num_err++;
		}
	}
	std::cout<<"\n number of error predicted is\t"<<num_err<<std::endl;
}

void learning_network::finalizing_multi_classes() {
	//relasing memory
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		delete[] experimental_results[n];
		delete[] targets_multi_class[n];
		delete[] outputsm[n];
		experimental_results[n] = nullptr;
		targets_multi_class[n] = nullptr;
		outputsm[n] = nullptr;
	}
	delete[] experimental_results;
	delete[] targets_multi_class;
	delete[] outputsm;

	std::size_t num_row = number_of_features * number_of_classes;
	for(std::size_t r = 0; r < num_row; r++) {
		delete[] weightsm[r];
		delete[] new_weightsm[r];
		delete[] gradient[r];
		delete[] delta_weight[r];
		delete[] hessian[r];
		delete[] Ihessian[r];
		weightsm[r] = nullptr;
		new_weightsm[r] = nullptr;
		gradient[r] = nullptr;
		delta_weight[r] = nullptr;
		hessian[r] = nullptr;
		Ihessian[r] = nullptr;
	}
	delete[] weightsm;
	delete[] new_weightsm;
	delete[] gradient;
	delete[] delta_weight;
	delete[] hessian;
	delete[] Ihessian;
	delete[] sum_w_experimental_results;
	delete[] predicted_output_multi_class;
	delete[] real_output;
	real_output = nullptr;
	predicted_output_multi_class = nullptr;
	sum_w_experimental_results = nullptr;
}
