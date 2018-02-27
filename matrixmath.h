#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _MATRIXMATH_H
#define _MATRIXMATH_H


typedef struct { double **m; int x, y;} matrix_m, *matrix;
typedef struct { matrix *A; double d; double *adj;} det_m, *detm;

extern matrix random_init(int x, int y, int range);
extern matrix matrix_cp(double *a, int x, int y);
extern matrix matrix_slice(matrix a, int x1, int x2, int y1, int y2);
extern void matrix_print(matrix a);
extern matrix cons_matrix(int x, int y, int number);
extern matrix eye_matrix(int x, int y, int number);
extern matrix matrix_t(matrix a);
extern matrix matrix_sum(matrix a, int dim);
extern matrix matrix_add(matrix a, matrix b);
extern matrix matrix_scal(matrix a, float f);
extern matrix matrix_sub(matrix a, matrix b);
extern matrix matrix_mul(matrix a, matrix b);
extern matrix matrix_miu(matrix a, int dim);
extern matrix matrix_cat(matrix a, matrix b, int dim);
extern double det(matrix A);
extern matrix matrix_inv(matrix a);
extern matrix matrix_reshape(matrix a, int x, int y);
extern matrix matrix_diag(matrix a);
extern matrix matrix_cov(matrix a);

#endif