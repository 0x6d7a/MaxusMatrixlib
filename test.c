#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixmath.h"
// #include "kpp_m.h"

void separate(void)
{
	printf("##########\n");
}

int main(int argc, char const *argv[])
{
	// matrix V = random_init(3, 1, 10);
	// matrix A = random_init(5, 3, 10);
	// point c;
	// c->v = A->m[0];
	// matrix B = cons_matrix(2, 2, 6);
	// matrix Br = random_init(2, 2, 10);
	// matrix C = eye_matrix(3, 3, 1);
	// matrix D = eye_matrix(3, 5, 2);
	// matrix sumA_col = matrix_sum(A, 1);
	// matrix sumA_row = matrix_sum(A, 2);

	double a[3][3] = {{1, 3, 2}, {-1, 0, 2}, {3, 1, -1}};
	matrix A = matrix_cp(&a[0][0], 3, 3);

	// matrix slice = matrix_slice(A, 2, 2, 1, 2);

	// matrix *AA = (matrix *)malloc(sizeof(matrix) * 3);
	// AA[0] = A;
	// AA[1] = B;
	// AA[2] = C;

	// matrix AplusC = matrix_add(A, C); // Done!
	// matrix As = matrix_scal(A, 2); // Done!
	// matrix AsubC = matrix_sub(A, C); // Done!
	// matrix BmBr = matrix_mul(B, Br); // Done!
	// matrix Amean_col = matrix_miu(A, 1); // Done!
	// matrix Amean_row = matrix_miu(A, 2); // Done!
	// matrix AT = matrix_t(A); // Done!
	// matrix AcatB = matrix_cat(A, V, 1); // Done!
	// matrix AcatB = matrix_cat(A, V, 2); // Done!
	// matrix Brespahe = matrix_reshape(B, 2, 2); // Done!
	// matrix Brespahe = matrix_reshape(V, 2, 2); // Done!
	// matrix Adiag = matrix_diag(A); // Done!
	// matrix Acov = matrix_cov(matrix_t(Aa));
	// matrix X = matrix_slice(A, 0, 2, 0, 0);

	// detm determinant = det(A, 1);
	//

	printf("Matrix A:\n");
	matrix_print(A);
	separate();

	matrix inv = matrix_inv(A);


	printf("Matrix Inv:\n");
	matrix_print(inv);
	separate();
	//
	// printf("%f\n", c->v[0]);
	// printf("%f\n", c->v[1]);
	// printf("%f\n", c->v[2]);

	// printf("Matrix A11:\n");
	// matrix_print(determinant->A[0]);
	// separate();

	// printf("Matrix A12:\n");
	// matrix_print(determinant->A[1]);
	// separate();

	// printf("Matrix A13:\n");
	// matrix_print(determinant->A[2]);
	// separate();

	// printf("det dA11:\n");
	// printf("%f\n", determinant->adj[0]);
	// separate();

	// printf("det dA12:\n");
	// printf("%f\n", determinant->adj[1]);
	// separate();

	// printf("det dA13:\n");
	// printf("%f\n", determinant->adj[2]);
	// separate();

	// printf("%f\n", determinant->d);


	// printf("Matrix AT:\n");
	// matrix_print(AT);
	// separate();

	// printf("Matrix B:\n");
	// matrix_print(B);
	// separate();

	// printf("Matrix C:\n");
	// matrix_print(C);
	// separate();

	// printf("Matrix Br:\n");
	// matrix_print(Br);
	// separate();

	// printf("Matrix D:\n");
	// matrix_print(D);
	// separate();

	// printf("Matrix sumA_col\n");
	// matrix_print(sumA_col);
	// separate();

	// printf("Matrix sumA_row\n");
	// matrix_print(sumA_row);
	// separate();

	// printf("Matrix Amean_col\n");
	// matrix_print(Amean_col);
	// separate();

	// printf("Matrix Acov\n");
	// matrix_print(Acov);
	// separate();

	return 0;
}
