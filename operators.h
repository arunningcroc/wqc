#ifndef OPERATORS
#define OPERATORS
#include <stdbool.h>
#include <complex.h>
#include "pool.h"
typedef struct Matrix {
    double complex **matrix;
    int rows;
    int cols;
} Matrix;
enum mat_types {
    EMPTY, CREATION, ANNIHILATION, IDENTITY
};
//hh
Matrix *tensor (pool *p, Matrix *a, Matrix *b);
Matrix *matmul (pool *p, Matrix *a, Matrix *b);
void print_matrix(Matrix *a);
bool cmpmat(Matrix *a, Matrix *b);
Matrix *matcreate(pool *p, int rows, int cols, enum mat_types mat_type);
Matrix *dagger (pool *p, Matrix *a);
Matrix **fock_basis (pool *p, int n);
complex double trace(Matrix *a);
double vnorm(Matrix *a, int col);
Matrix *matcoldiv(Matrix *a, int col, double div);
Matrix *matcolsub(Matrix *a, int col1, Matrix *b, int col2);
Matrix *matcolmul(Matrix *a, int col, double div);
Matrix *matcopy (pool *p, Matrix *a);
void matcopycol(Matrix *a, int col1, Matrix *b, int col2);
void qr (Matrix *a, Matrix *Q, Matrix *R);
Matrix **solve_eigenvalues(pool *p, Matrix *a);
Matrix *matsum(pool *p, Matrix *a, Matrix *b);
Matrix *matsub(pool *p, Matrix *a, Matrix *b);
Matrix *consmul(pool *p,Matrix *a, complex double c);
Matrix *trim_zero(pool *p, Matrix *a);
#endif

