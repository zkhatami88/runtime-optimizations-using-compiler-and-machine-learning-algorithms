#include <limits>
#include <math.h>
#include <iostream>

#define MAX_DOUBLE (std::numeric_limits<double>::max())

class learning_network {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	std::size_t number_of_classes;
	double** experimental_results; 				//the experimental values of the features of the training data
	double threshold; 							//the convergence for estimating the final weights

	//two-class logistic regression
	double* weightsb; 							//weights of our learning network
	double* new_weightsb;						//updated weights after each step	
	unsigned* targets_two_class;				//outputs of the training data
	double** diag_weightsb; 					//used for updating weights
	double* M; 									//used for updating weights
	
	//multi-class logistic regression
	double** weightsm; 							//weights of our learning network : (F * K) * 1
	double** new_weightsm;						//updated weights after each step : (F * K) * 1
	unsigned** targets_multi_class;				//real output of each experimental results : N * K
	double** outputsm;							//outputs of the training data : N * K
	double** gradient;							//gradient of E : (F * K) * 1
	double** hessian;							//hessian of E : (F * K) * (F * K)
	double** Ihessian;							//inverse of hessian of E : (F * K) * (F * K)
	double** delta_weight;						//used for updating weight : (F * K) * 1
	double sum_w_experimental_results;			//used for computing output


	
	void inverse_mat(double** mat, double** inverse, std::size_t size);
	
	/*
	two-class logistic regression (binary logistic regression) :
	*/
	void normalizing_weights_two_class();
	double get_value_of_M(std::size_t i);
	void updating_values_of_M_and_diag_weights();
	void get_cofactor(double **src, double **dest, std::size_t p, std::size_t q, std::size_t size);
	double calc_determinant( double **mat, std::size_t size);
	void adjoint(double** mat, double** adj, std::size_t size);	
	void inverse_martix_XSX(double** mat, double** inverse, std::size_t size);
	void updating_values_of_weights_two_class();
	void new_values_for_weightsb();
	double computing_new_least_squared_err_two_class();	
	void learning_weights_two_classes();
	void printing_weights_two_class();

	/*
	multi-class logistic regression 
	*/
	void normalizing_weights_multi_class();
	void convert_target_to_binary(unsigned* target_src, unsigned** targets_dst);
	void sum_w_experimental_results_n(std::size_t n);
	double output_nk(std::size_t n, std::size_t k);
	void computing_all_output_nk(std::size_t n, std::size_t k);
	void computing_all_output();
	int eye_kj(std::size_t k, std::size_t j);
	void gradient_E_k(std::size_t k);
	void computing_all_gradient();
	double hessian_E_kj(std::size_t k, std::size_t j);
	void computing_all_hessian();
	void computing_inverse_hessian();
	void learning_weights_multi_classes();
	void computing_delta_weight();
	void new_values_for_weightsm();
	double computing_new_least_squared_err_multi_class();	
	void updating_values_of_weights_multi_class();
	void printing_weights_multi_class();
	
public:
	learning_network(std::size_t number_of_expr, std::size_t number_of_ftrs, 
						std::size_t number_of_cls, double th) {
		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs;
		number_of_classes = number_of_cls;
		threshold = th;		
	}

	//initilializer for two_class
	void initilializer_two_class(double** expr_results, unsigned* target_expr){
		weightsb = new double[number_of_features +  1];
		new_weightsb = new double[number_of_features +  1];
		targets_two_class = new unsigned[number_of_experiments];
		diag_weightsb = new double*[number_of_experiments];
		M = new double[number_of_experiments];
		experimental_results = new double*[number_of_experiments];

		//initializing weights
		for(std::size_t i = 0; i < number_of_features +  1; i++) {
			weightsb[i] = 0.1;
		}				

		for(std::size_t i = 0; i < number_of_experiments; i++) {
			experimental_results[i] = new double[number_of_features + 1];
			targets_two_class[i] = target_expr[i];
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

	}

	//initilializer for multi_class
	void initilializer_multi_class(double** expr_results, unsigned* target_expr) {
		sum_w_experimental_results = 0;
		std::size_t num_row = number_of_features * number_of_classes;
		weightsm = new double*[num_row]; 							
		new_weightsm = new double*[num_row];
		gradient = new double*[num_row];
		delta_weight = 	new double*[num_row];
		hessian = new double*[num_row];
		Ihessian = new double*[num_row];
		experimental_results = new double*[number_of_experiments];
		targets_multi_class = new unsigned*[number_of_experiments];
		outputsm = new double*[number_of_experiments];

		for(std::size_t i = 0; i < num_row; i++) {
			weightsm[i] =  new double[1];
			new_weightsm[i] =  new double[1];
			gradient[i] =  new double[1];
			delta_weight[i] =  new double[1];
			hessian[i] =  new double[num_row];
			Ihessian[i] =  new double[num_row];

			//initializing weights
			weightsm[i][0] = 0.1;
		}

		
		for(std::size_t i = 0; i < number_of_experiments; i++) {
			experimental_results[i] = new double[number_of_features];
			targets_multi_class[i] = new unsigned[number_of_classes];
			outputsm[i] = new double[number_of_classes];

			//initializing experimental_results
			for(std::size_t f = 0; f < number_of_features; f++) {
				experimental_results[i][f] = expr_results[i][f];
			}
		}

		//initializing targets_multi_class
		convert_target_to_binary(target_expr, targets_multi_class);
	}

	//two-class logistic regression
	void learning_two_classes();
	double* retrieving_weights_two_classes();

	//multi-class logistic regression
	void learning_multi_classes();
	double** retrieving_weights_multi_classes();
};


/*
two-class logistic regression (binary logistic regression) :
*/

double learning_network::get_value_of_M(std::size_t i) {
	double temp = 0.0;
	for(std::size_t j = 0; j < number_of_features + 1; j++) {
		temp += weightsb[j] * experimental_results[i][j];
	}
	double result = (1 / (1 + double(1 / exp(temp))));
	return result;
}

void learning_network::updating_values_of_M_and_diag_weights() {
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		M[i] = get_value_of_M(i);
		diag_weightsb[i][i] = M[i] * (1 - M[i]);
	}
}

// calculate the cofactor of element (row,col)
void learning_network::get_cofactor(double **src, double **dest, std::size_t p, 
									std::size_t q, std::size_t size) {
	
	std::size_t i = 0, j = 0;

    // Looping for each element of the matrix
    for (std::size_t row = 0; row < size; row++) {
        for (std::size_t col = 0; col < size; col++) {
            //  Copying std::size_to temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q) {
                dest[i][j++] = src[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == size - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// Calculate the determinant recursively.
double learning_network::calc_determinant(double **mat, std::size_t size) {
	//size must be >= 0
	//stop the recursion when matrix is a single element
	double det = 0; // Initialize result
     if (size == 1) {
        return mat[0][0];
    }
 
 	// To store cofactors
    double** temp = new double*[size];
    for(std::size_t i = 0; i < size; i++) {
    	temp[i] = new double[size];
    } 
 
    double sign = 1.0;  // To store sign multiplier
 
     // Iterate for each element of first row
    for (std::size_t f = 0; f < size; f++) {
        // Getting Cofactor of mat[0][f]
        get_cofactor(mat, temp, 0, f, size);
        det += sign * mat[0][f] * calc_determinant(temp, size - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    } 

    for(std::size_t i = 0; i < size; i++) {
    	delete[] temp[i];
    	temp[i] = nullptr;
    }

    delete[] temp;

    return det;
}

// Function to get adjoint of mat in adj
void learning_network::adjoint(double** mat, double** adj, std::size_t size) {
    if (size == 1) {
        adj[0][0] = 1;
        return;
    }
 
    // temp is used to store cofactors of mat
    double sign = 1.0;
    double** temp = new double*[size];
    for(std::size_t i = 0; i < size; i++) {
    	temp[i] = new double[size];
    }
 
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            // Get cofactor of mat[i][j]
            get_cofactor(mat, temp, i, j, size);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1: -1;
 
            // std::size_terchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(calc_determinant(temp, size - 1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
void learning_network::inverse_mat(double** mat, double** inverse, std::size_t size) {    
    double det = calc_determinant(mat, size);
    std::cout<<"\n det is : "<<det<<std::endl;
    if (det == 0) {
        std::cout << "Singular matrix, can't find its inverse";
        return;
    }
 
    // Find adjoint
    double** adj = new double*[size];
    for(std::size_t i = 0; i < size; i++) {
    	adj[i] = new double[size];
    }

    adjoint(mat, adj, size);
 
    // Find Inverse using formula "inverse(mat) = adj(mat)/det(mat)"
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            inverse[i][j] = adj[i][j]/det;
        }
    }

	for(std::size_t i = 0; i < size; i++) {
		delete adj[i];
	}    
	delete[] adj;
 }

void learning_network::inverse_martix_XSX(double** mat, double** inverse, std::size_t size) {
    inverse_mat(mat, inverse, size);
}

void learning_network::new_values_for_weightsb() {
	//allocating memory
	double** temp1; //XSX : f * f
	double** temp2; //SXW + y - M : expr * 1
	double** temp3; //XSX * X^T : f * expr
	double** temp4; //X^T * S
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
			//initializing value
			if(j == 0) {
				temp4[i][j] = 0;	
			}
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
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		delete temp1[i];
		delete temp3[i];
		temp1[i] = nullptr;
		temp3[i] = nullptr;
	}
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		delete temp2[i];
		temp2[i] = nullptr;
	}
	delete[] temp1;
	delete[] temp2;
	delete[] temp3;
}

double learning_network::computing_new_least_squared_err_two_class() {
	double temp_err = 0.0;
	double old_wx, new_wz;
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		old_wx = 0;
		new_wz = 0;
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			old_wx += weightsb[j] * experimental_results[i][j];
			new_wz += new_weightsb[j] * experimental_results[i][j];
		}
		temp_err += pow( ((1 / (1 + (1 / exp(old_wx))))) - (1 / (1 + (1 / exp(new_wz)))) , 2.0);
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
	double least_squared_err = MAX_DOUBLE;	
	std::size_t itr = 0;
	while(threshold < least_squared_err) {		
		updating_values_of_M_and_diag_weights();
		new_values_for_weightsb();
		least_squared_err = computing_new_least_squared_err_two_class();
		std::cout<<"("<<itr<<")\t"<<"Least_squared_err =\t" << least_squared_err<<std::endl;
		updating_values_of_weights_two_class();		
		printing_weights_two_class();
		itr++;
	}
}

//FIXME
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

	std::cout<<"The values of training data are : "<<std::endl;
	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			std::cout<<experimental_results[i][j]<<" , ";
		}
		std::cout<<std::endl;
	}


	//releasing memory
	delete[] averages;
	delete[] averages_2;
	delete[] var;
}

void learning_network::learning_two_classes() {
	//normalizing_weights_two_class();
	learning_weights_two_classes();
}

double* learning_network::retrieving_weights_two_classes() {
	return weightsb;
}

/*
multi-class logistic regression 
*/

void learning_network::convert_target_to_binary(unsigned* target_src, unsigned** targets_dst) {
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
int learning_network::eye_kj(std::size_t k, std::size_t j) {
	if(k == j) {
		return 1;
	}
	return 0;
}

//computing ouput
void learning_network::sum_w_experimental_results_n(std::size_t n) {	
	for(std::size_t k = 0; k < number_of_classes; k++) {
		std::size_t offset = k * number_of_features;
		for(std::size_t f = 0; f < number_of_features; f++) {
			sum_w_experimental_results += weightsm[offset + f][0] * experimental_results[n][f];
		}
	}
}

double learning_network::output_nk(std::size_t n, std::size_t k) {
	double temp = 0.0;
	std::size_t offset = k * number_of_features;
	for(std::size_t f = 0; f < number_of_features; f++) {
		temp += weightsm[offset + f][0] * experimental_results[n][f];
	}

	double out = exp(temp)/exp(sum_w_experimental_results); 
	return out;
}

void learning_network::computing_all_output_nk(std::size_t n, std::size_t k) {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		sum_w_experimental_results = 0;
		sum_w_experimental_results_n(n);
		for(std::size_t k = 0; k < number_of_classes; k++) {
			outputsm[n][k] = output_nk(n, k);
			//std::cout<<outputsm[n][k]<<"\t";
		}
		//std::cout<<std::endl;
	}
}

void learning_network::computing_all_output() {
	for(std::size_t n = 0; n < number_of_experiments; n++) {
		for(std::size_t k = 0; k < number_of_classes; k++) {
			computing_all_output_nk(n, k);
		}
	}
}

//computing gradient for each W
void learning_network::gradient_E_k(std::size_t k) {
	std::size_t offset = k * number_of_features;
	for(std::size_t f = 0; f < number_of_features; f++) {
		gradient[offset + f][0] = 0;
		for(std::size_t n = 0; n < number_of_experiments; n++) {
			gradient[offset + f][0] += (output_nk(n, k) - targets_multi_class[n][k]) * experimental_results[n][f];
		}
	}
}

void learning_network::computing_all_gradient(){
	for(std::size_t k = 0; k < number_of_classes; k++) {
		gradient_E_k(k);
	}
	for(std::size_t i = 0; i < number_of_classes * number_of_features; i++) {
		//std::cout<<gradient[i][0]<<",";
	}
	//std::cout<<std::endl;
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

	std::size_t num_row = number_of_features * number_of_classes;
	for(std::size_t i = 0; i < num_row; i++) {
		for(std::size_t j = 0; j < num_row; j++) {
			//std::cout<<hessian[i][j]<<"\t";
		}
		//std::cout<<std::endl;
	}
}

void learning_network::computing_inverse_hessian() {
	std::size_t size = number_of_classes * number_of_features;
	inverse_mat(hessian, Ihessian, size);
	//writing output of invesring hessian in text file
	for(std::size_t i = 0; i < size; i++) {
		for(std::size_t j = 0; j < size; j++) {
			std::cout<<Ihessian[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
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
		new_weightsm[i][0] = weightsm[i][0];
	}
}

void learning_network::printing_weights_multi_class() {
	std::size_t size = number_of_classes * number_of_features;

	for(std::size_t k = 0; k < number_of_classes; k++) {
		for(std::size_t f = 0; f < number_of_features; f++) {
			std::size_t offset = k * number_of_features + f;
			std::cout<<"weights["<<k<<"]["<<f<<"] = "<<weightsm[offset][0]<<"\t";
		}
		std::cout<<std::endl;
	}
	std::cout<<"\n----------\n"<<std::endl;
}

//updating weights till error meets the defined threshold
void learning_network::learning_weights_multi_classes() {
	double least_squared_err = MAX_DOUBLE;

	//while(threshold < least_squared_err) {
		computing_all_output();		
		computing_all_gradient();
		computing_all_hessian();		
		new_values_for_weightsm();
		/*
		least_squared_err = computing_new_least_squared_err_multi_class();
		std::cout<<"\nLeast_squared_err:" << least_squared_err<<std::endl;
		updating_values_of_weights_multi_class();
		printing_weights_multi_class();*/
	//}
}

void learning_network::normalizing_weights_multi_class() {
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

	std::cout<<"The values of training data are : "<<std::endl;
	//computing average and variance values for each feature
	for(std::size_t i = 0; i < number_of_experiments; i++) {		
		for(std::size_t j = 0; j < number_of_features; j++) {
			std::cout<<experimental_results[i][j]<<" , ";
		}
		std::cout<<std::endl;
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
