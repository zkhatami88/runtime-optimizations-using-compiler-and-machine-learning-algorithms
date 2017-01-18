#include <limits>
#include <math.h>
#include <iostream>

#define MAX_DOUBLE (std::numeric_limits<double>::max())

class learning_network {

	std::size_t number_of_experiments;
	std::size_t number_of_features;
	double* weights; 							//weights of our learning network
	double* new_weights;						//updated weights after each step
	double** experimental_results; 				//the experimental values of the features of the training data
	double* outputs;							//outputs of the training data
	double** diag_weights; 						//used for updating weights
	double* M; 									//used for updating weights
	double threshold; 							//the convergence for estimating the final weights

	
	double get_value_of_M(std::size_t i);
	void updating_values_of_M_and_diag_weights();
	void get_cofactor(double **src, double **dest, std::size_t p, std::size_t q, std::size_t size);
	double calc_determinant( double **mat, std::size_t size);
	void adjoint(double** mat, double** adj, std::size_t size);
	void inverse_mat(double** mat, double** inverse, std::size_t size);
	void inverse_martix_XSX(double** mat, double** inverse, std::size_t size);
	void updating_values_of_weights();
	void new_values_for_weights();
	double computing_new_least_squared_err();
	void normalizing_weights();
	void learning_weights();	
	
public:
	learning_network(std::size_t number_of_expr, std::size_t number_of_ftrs, 
						double** expr_results, double* output_expr, double th) {

		number_of_experiments = number_of_expr;
		number_of_features = number_of_ftrs; 

		weights = new double[number_of_features +  1];
		new_weights = new double[number_of_features +  1];
		for(std::size_t i = 0; i < number_of_features +  1; i++) {
			weights[i] = 0.1;
		}

		experimental_results = new double*[number_of_experiments];
		outputs = new double[number_of_experiments];
		diag_weights = new double*[number_of_experiments];
		M = new double[number_of_experiments];

		for(std::size_t i = 0; i < number_of_experiments; i++) {
			experimental_results[i] = new double[number_of_features + 1];
			outputs[i] = output_expr[i];
			diag_weights[i] = new double[number_of_experiments];
			M[i] = 0.0;

			experimental_results[i][0] = 1.0;
			for(std::size_t j = 0; j < number_of_features + 1; j++) {
				experimental_results[i][j] = expr_results[i][j];
			}
			for(std::size_t j = 0; j < number_of_experiments; j++) {
				diag_weights[i][j] = 0.0;
			}
		}

		threshold = th;
	}

	void learning();
	double* retrieving_weights();
};

double learning_network::get_value_of_M(std::size_t i) {
	double temp = 0.0;
	for(std::size_t j = 0; j < number_of_features + 1; j++) {
		temp += weights[j] * experimental_results[i][j];
	}
	double result = (1 / (1 + double(1 / exp(temp))));
	return result;
}

void learning_network::updating_values_of_M_and_diag_weights() {
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		M[i] = get_value_of_M(i);
		diag_weights[i][i] = M[i] * (1 - M[i]);
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

void learning_network::new_values_for_weights() {
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
			temp4[i][j] = experimental_results[j][i] * diag_weights[j][j];
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
			temp2[i][0] += diag_weights[i][i] * experimental_results[i][j] * weights[j]; 
		}
		temp2[i][0] += outputs[i] - M[i];
	}

	//updating weights
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		new_weights[i] = 0.0;
		for(std::size_t j = 0; j < number_of_experiments; j++) {
			for(std::size_t k = 0; k < number_of_features + 1; k++) {
				if(k == 0) {
					temp3[i][j] = 0.0;
				}
				temp3[i][j] += inverse_temp1[i][k] * experimental_results[j][k];
			}
			new_weights[i] += temp3[i][j] * temp2[j][0];
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

double learning_network::computing_new_least_squared_err() {
	double temp_err = 0.0;
	double old_wx, new_wz;
	for(std::size_t i = 0; i < number_of_experiments; i++) {
		old_wx = 0;
		new_wz = 0;
		for(std::size_t j = 0; j < number_of_features + 1; j++) {
			old_wx += weights[j] * experimental_results[i][j];
			new_wz += new_weights[j] * experimental_results[i][j];
		}
		temp_err += pow( ((1 / (1 + (1 / exp(old_wx))))) - (1 / (1 + (1 / exp(new_wz)))) , 2.0);
	}
	return sqrt(temp_err);
}

void learning_network::updating_values_of_weights() {
	for(std::size_t i = 0; i < number_of_features + 1; i++) {
		weights[i] = new_weights[i];
	}
}

//this func applys experimental values of the training data on our learning network
void learning_network::learning_weights() {
	double least_squared_err = MAX_DOUBLE;	

	while(threshold < least_squared_err) {	
		std::cout <<"\n ====================\n"<<std::endl;	
		updating_values_of_M_and_diag_weights();
		new_values_for_weights();
		least_squared_err = computing_new_least_squared_err();
		std::cout<<"\n least_squared_err:" << least_squared_err<<std::endl;
		updating_values_of_weights();		
		for(std::size_t i = 0; i < number_of_features + 1; i++) {
			std::cout<<"weights["<<i<<"] = "<<weights[i]<<" , ";
		}
		std::cout<<std::endl;
	}
}

void learning_network::normalizing_weights() {
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
}

void learning_network::learning() {
	//normalizing_weights();
	learning_weights();
}

double* learning_network::retrieving_weights() {
	return weights;
}
