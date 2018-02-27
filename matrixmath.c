#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixmath.h"

matrix random_init(int x, int y, int range)
{
	int i, j;

#ifndef for_x
#define for_x for(i; i < x; i++)
#endif

#ifndef for_y
#define for_y for(j; j < y; j++)
#endif

	i = 0;
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * y);
	}

	result->x = x;
	result->y = y;

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = rand() % range;
		}
	}

	return result;
	free(result);
}

matrix matrix_cp(double *a, int x, int y)
{
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	result->x = x;
	result->y = y;

	int i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * y);
		result->m[i] = a + i * y;
	}

	return result;
	free(result);
}

matrix matrix_slice(matrix a, int x1, int x2, int y1, int y2)
{
	/**
	 * Parameters: 
	 * x1 - starting row
	 * x2 - ending (including) row
	 * y1 - starting col
	 * y2 - ending (including) col
	 */
	
	int i, j;
	int x, y;
	x = x2 - x1 + 1;
	y = y2 - y1 + 1;

	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->x = x;
	result->y = y;
	result->m = (double **)malloc(sizeof(double *) * x);

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * y);
		result->m[i] = a->m[x1 + i] + y1;
	}
	return result;
	free(result);
}

void matrix_print(matrix a)
{
	int x, y;
	int i, j;

	x = a->x;
	y = a->y;

	printf("Size [%d,%d]\n", x, y);

	char next;

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			if( j != y - 1)
			{
				next = '\t';
			}
			else
			{
				next = '\n'; 
			}
			printf("%f%c", a->m[i][j], next);
		}
	}
}

matrix cons_matrix(int x, int y, double number)
{

	int i, j;

	i = 0;
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * y);
	}

	result->x = x;
	result->y = y;

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = number;
		}
	}

	return result;
	free(result);
}

matrix eye_matrix(int x, int y, int number)
{
	int i, j;

	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * y);
	}

	result->x = x;
	result->y = y;

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			if(i != j)
			{
				result->m[i][j] = 0;
				continue;	
			}
			result->m[i][j] = number;
		}
	}

	return result;
	free(result);
}

matrix matrix_t(matrix a)
{

	int x = a->x;
	int y = a->y;
	int i, j;

	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * y);

	result->x = y;
	result->y = x;

	j = 0;
	for_y
	{
		result->m[j] = (double *)malloc(sizeof(double) * x);
	}

	j = 0;
	for_y
	{
		i = 0;
		for_x
		{
			result->m[j][i] = a->m[i][j];
		}
	}

	return result;
	free(result);
}

matrix matrix_sum(matrix a, int dim)
{
	int x = a->x;
	int y = a->y;
	int i, j;
	double sum;

	matrix result;

	if(dim == 1)
	{
		result = (matrix_m *)malloc(sizeof(matrix_m));
		result->x = x;
		result->y = 1;
		i = 0;
		
		result->m = (double **)malloc(sizeof(double *) * x);
		for_x
		{
			sum = 0;
			result->m[i] = (double *)malloc(sizeof(double));

			j = 0;
			for_y
			{
				sum += a->m[i][j];
			}
			result->m[i][0] = sum;
		}
	}
	else if(dim == 2)
	{
		matrix t = matrix_t(a);
		result = matrix_t(matrix_sum(t, 1));
	}

	return result;
	free(result);
}

matrix matrix_add(matrix a, matrix b)
{

	if((a->x != b->x) || (a->y != b->y))
	{
		printf("Dimensionality mismatch.\n");
		exit(0);
	}

	int x = a->x;
	int y = a->y;
	int i, j;
	double sum;
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	result->x = x;
	result->y = y;

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) *y);
	}

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = a->m[i][j] + b->m[i][j];
		}
	}

	return result;
	free(result);
}

matrix matrix_scal(matrix a, float f)
{
	int x = a->x;
	int y = a->y;
	int i, j;
	
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	result->x = x;
	result->y = y;

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) *y);
	}

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = a->m[i][j] * f;
		}
	}

	return result;	
	free(result);
}

matrix matrix_sub(matrix a, matrix b)
{
	matrix c = matrix_scal(b, -1);

	return matrix_add(a, c);
}

matrix matrix_mul(matrix a, matrix b)
{

	if(a->y != b->x)
	{
		printf("Dimensionality mismatch.");
		exit(0);
	}

	int x = a->x;
	int y = a->y;
	int z = b->y;
	int i, j, k;

	double sum;

#ifndef for_z
#define for_z for(k; k < z; k++)
#endif

	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *) * x);

	result->x = x;
	result->y = z;

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * z);
	}

	i = 0;
	for_x
	{
		k = 0;
		for_z
		{
			j = 0;
			sum = 0;

			for_y
			{
				sum += a->m[i][j] * b->m[j][k];
			}
			result->m[i][k] = sum;
		}
	}

	return result;
	free(result);
}

matrix matrix_miu(matrix a, int dim)
{
	int x = a->x;
	int y = a->y;
	int i, j;
	matrix sum;

	if(dim == 1)
	{
		sum = matrix_sum(a, 1);
		i = 0;
		for_x
		{
			sum->m[i][0] /= y;
		}

	}
	else if(dim == 2)
	{
		sum = matrix_sum(a, 2);
		j = 0;
		for_y
		{
			sum->m[0][j] /= x;
		}
	}

	return sum;
}

matrix matrix_cat(matrix a, matrix b, int dim)
{
	int x, y;
	int tx, ty;
	int i, j;

	if(dim == 1)
	{
		tx = a->x + b->x;
		ty = fmax(a->y, b->y);
	}
	else if(dim == 2)
	{
		tx = fmax(a->x, b->x);
		ty = a->y + b->y;
	}
	else
	{
		printf("Dimensionality error."); 
		exit(0);
	}

	matrix result = cons_matrix(tx, ty, 0);

	i = 0;
	x = a->x;
	y = a->y;
	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = a->m[i][j];
		}
	}

	x = b->x;
	y = b->y;
	i = 0;
	for_x
	{
		j = 0;
		for_y
		{	
			if(dim == 1) result->m[a->x+i][j] = b->m[i][j];
			else if(dim == 2) result->m[+i][a->y + j] = b->m[i][j];
		}
	}

	return result;
	free(result);

}

double det(matrix a)
{
	int x = a->x, y = a->y;

	if( x != y)
	{
		printf("Dimensionality error.\n");
		exit(0);
	}
	else if(x == 1)
	{
		return a->m[0][0];
	}

	double sum = 0;
	matrix aij;
	int i, j;
	i = 0;
	j = 0;

	for_y
	{
		if(j == 0)
		{
			aij = matrix_slice(a, i + 1, x - 1, j + 1, y - 1);
			sum += a->m[i][j] * pow(-1, i + j) * det(aij);
		}
		else if( j < y - 1)
		{
			matrix ajl = matrix_slice(a, i + 1, x - 1, 0, j - 1);
			matrix ajr = matrix_slice(a, i + 1, x - 1, j + 1, y - 1);
			aij = matrix_cat(ajl, ajr, 2);
			sum += a->m[i][j] * pow(-1, i + j) * det(aij);
		}
		else
		{
			aij = matrix_slice(a, i + 1, x - 1, 0, j - 1);
			sum += a->m[i][j] * pow(-1, i + j) * det(aij);
		}
	}
	return sum;
}

matrix matrix_inv(matrix a)
{
	int i, j, x, y;
	i = 0;

	if(det(a) == 0.0)
	{
		printf("Uninversable matrix | Determinant is 0.\n");
		exit(0);
	}

	if(a->x == 1)
	{
		return matrix_scal(a, 1.0/det(a));
	}

	x = a->x;
	y = a->y;
	matrix adjoint = cons_matrix(a->x, a->y, 0);
	matrix aij, ajl, ajr, asrow_a, asrow_b, ascol_a, ascol_b, ar;

	for_x
	{
		if(i == 0)
		{
			j = 0;
			for_y
			{
				if(j == 0)
				{
					aij = matrix_slice(a, i + 1, x - 1, j + 1, y - 1);
				}
				else if( j < y - 1)
				{
					ajl = matrix_slice(a, i + 1, x - 1, 0, j - 1);
					ajr = matrix_slice(a, i + 1, x - 1, j + 1, y - 1);
					aij = matrix_cat(ajl, ajr, 2);
				}
				else
				{
					aij = matrix_slice(a, i + 1, x - 1, 0, j - 1);
				}
				adjoint->m[i][j] = pow(-1, i + j) * det(aij);
			}
		}
		else if( i < x - 1)
		{
			j = 0;
			for_y
			{
				asrow_a = matrix_slice(a, 0, i - 1, 0, y - 1);
				asrow_b = matrix_slice(a, i + 1, x - 1, 0, y - 1);
				ar = matrix_cat(asrow_a, asrow_b, 1);

				if(j == 0)
				{
					aij = matrix_slice(ar, 0, ar->x - 1, j + 1, y - 1);
				}
				else if( j < y - 1)
				{
					ajl = matrix_slice(ar, 0, ar->x - 1, 0, j - 1);
					ajr = matrix_slice(ar, 0, ar->x - 1, j + 1, y - 1);
					aij = matrix_cat(ajl, ajr, 2);
				}
				else
				{
					aij = matrix_slice(ar, 0, ar->x - 1, 0, j - 1);
				}
				adjoint->m[i][j] = pow(-1, i + j) * det(aij);
			}
		}
		else
		{
			j = 0;
			for_y
			{
				if(j == 0)
				{
					aij = matrix_slice(a, 0, i - 1, j + 1, y - 1);
				}
				else if( j < y - 1)
				{
					ajl = matrix_slice(a, 0, i - 1, 0, j - 1);
					ajr = matrix_slice(a, 0, i - 1, j + 1, y - 1);
					aij = matrix_cat(ajl, ajr, 2);
				}
				else
				{
					aij = matrix_slice(a, 0, i - 1, 0, j - 1);
				}

				adjoint->m[i][j] = pow(-1, i + j) * det(aij);
			}
		}
	}

	matrix result = matrix_scal(matrix_t(adjoint), 1.0/det(a));

	return result;

	free(result);
	free(ajl);
	free(ajr);
	free(aij);
	free(asrow_a);
	free(asrow_b);
	free(ascol_a);
	free(ascol_b);
	free(ar);
}

matrix matrix_reshape(matrix a, int x, int y)
{
	int i, j;

	matrix result_col = a;

	i = 1;
	for_x
	{
		result_col = matrix_cat(result_col, a, 1);
	}

	matrix result = result_col;

	j = 1;
	for_y
	{
		result = matrix_cat(result, result_col, 2);
	}

	result->x = x * a->x;
	result->y = y * a->y;

	return result;
	free(result);
	free(result_col);
}

matrix matrix_diag(matrix a)
{
	int x = a->x;
	int y = a->y;
	int i, j;
	
	matrix result = (matrix_m *)malloc(sizeof(matrix_m));
	result->m = (double **)malloc(sizeof(double *));

	result->x = 1;
	result->y = fmin(x, y);

	i = 0;
	for_x
	{
		result->m[i] = (double *)malloc(sizeof(double) * fmin(x, y));
	}

	i = 0;
	for_x
	{
		j = 0;
		for_y
		{
			if(i != j ) continue;
			result->m[0][j] = a->m[i][j];
		}
	}

	return result;
	free(result);
}


matrix matrix_cov(matrix a)
{
	matrix miu = matrix_reshape(matrix_miu(a, 1), 1, a->y);
	matrix decentered = matrix_sub(a, miu);
	matrix result = matrix_mul(decentered, matrix_t(decentered));
	
	return matrix_scal(result, (float)1.0/(a->y - 1));
	free(result);
	free(miu);
	free(decentered);
}


matrix matrix_log(matrix a)
{
	matrix result = cons_matrix(a->x, a->y, 0);

	int i, j, x, y;
	i = 0;
	x = a->x;
	y = a->y;

	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = log(a->m[i][j]);
		}
	}

	return result;
	free(result);
}

matrix matrix_square(matrix a)
{
	matrix result = cons_matrix(a->x, a->y, 0);

	int i, j, x, y;
	i = 0;
	x = a->x;
	y = a->y;

	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = pow(a->m[i][j], 2);
		}
	}

	return result;
	free(result);
}

matrix matrix_exp(matrix a)
{
	matrix result = cons_matrix(a->x, a->y, 0);

	int i, j, x, y;
	i = 0;
	x = a->x;
	y = a->y;

	for_x
	{
		j = 0;
		for_y
		{
			result->m[i][j] = exp(a->m[i][j]);
		}
	}

	return result;
	free(result);
}

matrix matrix_x(matrix a, matrix b)
{
	int i, j, x, y;

	if(a->x == b->x && b->y == 1)
	{
		matrix_x(a, matrix_reshape(b, 1, a->y));
	}
	else if(a->x == b->x && a->y == b->y)
	{
		matrix result = cons_matrix(a->x, a->y, 0);
		i = 0;
		x = a->x;
		y = a->y;

		for_x
		{
			j = 0;
			for_y
			{
				result->m[i][j] = a->m[i][j] * b->m[i][j];
			}
		}

		return result;
		free(result);
	}
	else if(a->x == 1 && a->y == b->y)
	{
		matrix_x(matrix_reshape(a, b->x, 1), b);
	}
	else
	{
		printf("Dimensionality error.\n");
		exit(0);
	}

}


matrix matrix_d(matrix a, matrix b)
{
	int i, j, x, y;

	if(a->x == b->x && a->y == b->y)
	{
		matrix result = cons_matrix(a->x, a->y, 0);
		i = 0;
		x = a->x;
		y = a->y;

		for_x
		{
			j = 0;
			for_y
			{
				result->m[i][j] = (double) a->m[i][j] / b->m[i][j];
			}
		}

		return result;
		free(result);
	}
	else
	{
		printf("Dimensionality error.\n");
		exit(0);
	}

}