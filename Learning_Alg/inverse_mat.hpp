
//calculate transpose of matrix
void transpose(double** c, double** dst, int size, double determinte) {	
	double** b = new double*[size];
	for(unsigned i = 0; i < size; i++) {
		b[i] = new double[size];
	}
	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = 0; j < size; j++) {
			b[i][j] = c[j][i];
		}
	}
	for (unsigned i = 0; i < size; i++){
		for (unsigned j = 0; j < size; j++){
			dst[i][j] = b[i][j]/determinte;
		}
	}
}

//calculate cofactor of matrix
void get_cofactor(double** mat, double** dst ,int size, double determinte) {
	double** b = new double*[size];
	double** c = new double*[size];
	for(unsigned i = 0; i < size; i++) {
		b[i] = new double[size];
		c[i] = new double[size];
	}
	for (unsigned h = 0; h < size; h++) {
		for (unsigned l = 0;l < size; l++) {
			unsigned m=0, k=0;
			for (unsigned i = 0; i < size; i++) {
				for (unsigned j = 0; j < size; j++) {
					if (i != h && j != l){
						b[m][k] = mat[i][j];
						if (k < (size - 2)){
							k++;
						}
						else{
							k = 0;
							m++;
						}
					}
				}
			}
			c[h][l] = pow(-1, (h + l)) * calc_determinant(b, (size - 1));
		}
	}
	transpose(c, dst, size, determinte);
}

//calculate minor of matrix OR build new matrix : k-had = minor
void minor(double** b, double** a, unsigned i, int size) {
	unsigned h = 0, k = 0;
	for(unsigned l = 1; l < size; l++) {
		for(unsigned j = 0; j < size; j++){
			if(j == i){
				continue;
			}
			b[h][k] = a[l][j];
			k++;
			if(k == (size - 1)){
				h++;
				k=0;
			}
		}
	}
}

//calculate determinte of matrix
double calc_determinant(double** mat, int size) {
	double** temp = new double*[size];
	for(unsigned i = 0; i < size; i++) {
		temp[i] = new double[size];
	}
	double determinte = 0.0;
	if (size == 1) {
		return mat[0][0];
	}
	else if(size == 2) {
		return (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);
	}
	else {
		for(unsigned i=0; i < size; i++){
			minor(temp, mat, i, size);
			determinte = (double) (determinte + mat[0][i] * pow(-1,i) * det(temp, (size-1)));
		}
	}

	return determinte;
}

void inverse_mat(double** mat, double** inverse_mat, int size, double determinte) {
	if(determinte == 0) {
		printf("\nThis mat is not invertable\n");
	}
	else if(size == 1) {
		inverse_mat[0][0] = 1;
	}
	else {
		get_cofactor(mat, inverse_mat, size, determinte);
	}
}