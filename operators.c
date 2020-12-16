#include "operators.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
/*Tensor of two operators living in different Hilbert spaces - so, a Kronecker product.*/
Matrix *
tensor(pool   *p,
       Matrix *a,
       Matrix *b)
{

    Matrix *res = matcreate(p, a->rows*b->rows, a->cols*b->cols, EMPTY);
    for(int ai = 0; ai < a->rows; ai++) {
        for(int aj = 0; aj < a->cols; aj++) {
            for(int i = 0; i < b->rows; i++) {
                for(int j = 0; j < b->cols; j++) {
                    res->matrix[ai*b->rows + i][aj*b->cols + j] = a->matrix[ai][aj] * b->matrix[i][j];
                }
            }
        }
    }
    return res;
}
/*Calculates the trace of a square matrix*/
complex double
trace(Matrix *a)
{   
    complex double d = 0.0;
    if (a->rows != a->cols) printf("Incorrect matrix dimensions for trace");
    for(int i = 0; i < a->rows; i++) {
        d += a->matrix[i][i];
    }
    return d;
}
/* fock_basis generates a number basis of given dimension n, with |0> = (1,0,0,0..)^T */
Matrix **
fock_basis(pool *p,
           int n)
{
    Matrix **a = pool_alloc(p, n * sizeof(Matrix *));
    for(int i = 0; i < n; i++) {
        a[i] = matcreate(p, n, 1, EMPTY);
        a[i]->matrix[i][0] = 1;
    }
    return a;
}
/* matmul is a standard matrix multiplication routine.
 * there's nothing special about it, except perhaps that
 * it is very slow. matmul_inplace multiplies wihtout
 * allocating more memory - result is stored on the left hand side.
 * that's an optimization trick for eigenvalue solving.
 */

Matrix *
matmul (pool *p,
         Matrix *a,
         Matrix *b)
{
    if(a->cols != b->rows) {
        printf("Wrong matrix dimensions");
        abort();
    }
    complex double sum = 0.0;
    enum mat_types t = EMPTY;
    Matrix *prod = matcreate(p, a->rows, b->cols, t);
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < b->cols; j++) {
            for(int k = 0; k < b->rows; k++) {
                sum = sum + a->matrix[i][k]*b->matrix[k][j];
            }
            prod->matrix[i][j] = sum;
            sum = 0;
        }
    }
    return prod;
}
/*Calculates the QR decomposition of a matrix.*/
void
qr (Matrix *a,
    Matrix *Q,
    Matrix *R)
{

    pool *p = begin_calculation(10000*sizeof(double));
    Matrix *T = matcreate(p, a->rows, 1, EMPTY);
    Matrix *S = matcreate(p, a->rows, 1, EMPTY);
    for (int i = 0; i < a->cols; i++) {
        //print_matrix(Q);
        matcopycol(a,i,Q,i);
        for (int j = 0; j < i; j++) {
            //printf("%f\n", creal(T->matrix[0][0]));
            matcopycol(Q,j,T,0);
            //printf("%f\n", creal(T->matrix[0][0]));
            matcopycol(a,i,S,0);
            
            complex double r = 0;
            for(int k = 0; k < a->rows; k++) {
                r += T->matrix[k][0] * S->matrix[k][0];
            }
            //printf("%f\n", creal(r));
            R->matrix[j][i] = r;
            matcolsub(Q, i, matcolmul(T,0,r), 0);
            
        } 
        R->matrix[i][i] = vnorm(Q,i);
        matcoldiv(Q,i,R->matrix[i][i]);
    }
    end_calculation(p);
}
/*Lenght of a column vector in the matrix.*/
double vnorm(Matrix *a,
             int col)
             
{   
    double sum = 0.0;
    for(int i = 0; i < a->rows; i++) sum += a->matrix[i][col]*a->matrix[i][col];
    return sqrt(sum);
    
}
Matrix  **
solve_eigenvalues(pool   *p,
                  Matrix *a)
                  
{
    int counter = 0;
    int iteration = 0;
    double target = 1e-2;
    Matrix *V = matcreate(p, a->rows, a->cols, IDENTITY);
    Matrix *Q = matcreate(p, a->rows, a->cols, EMPTY);
    Matrix *R = matcreate(p, a->rows, a->cols, EMPTY);
    Matrix **eigs = pool_alloc(p, 2 * sizeof(Matrix *));
    //print_matrix(a);
    while(1) {
        qr(a, Q, R);
        a = matmul(p, R, Q);
        V = matmul(p, V, Q);
        for(int i = 0; i < a->rows ; i++) {
            for(int j = 0; j < a->cols; j++) {
                if (i != j && fabs(creal(a->matrix[i][j])) < target && fabs(cimag(a->matrix[i][j])) < target)  {
                    ++counter;
                }
            }
        }
        //print_matrix(a);
        if(counter == a->rows * a->cols - a->rows) {
            printf("Found eigenvalues");
            break;
        }
        else if(iteration > 5) {
            printf("Eigensolver did not converge for some reason. You're shit outta luck, boy. Here's the final matrix.\n");
            //print_matrix(a);
            break;
        }
        else {
            counter = 0;
            iteration++;
        }
    }
    /*for(int i = 0; i < a->rows; i++) {
        eigenvals[i] = a->matrix[i][i];
    }*/
    eigs[0] = V;
    eigs[1] = a;
    return eigs;
    
}
/*Sum the two matrices a and b*/
Matrix *
matsum(pool   *p,
       Matrix *a,
       Matrix *b)
{
    if(a->cols != b->cols || a->rows != b->rows) {
        printf("Dimension error in matsum");
        return NULL;
    }
    Matrix *c = matcreate(p, a->rows, a->cols, EMPTY);
    
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            c->matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];
        }
    }
    return c;
}
//subtract b from a.
Matrix *
matsub(pool   *p,
       Matrix *a,
       Matrix *b)
{
    if(a->cols != b->cols || a->rows != b->rows) {
        printf("Dimension error in matsum");
        return NULL;
    }
    Matrix *c = matcreate(p, a->rows, a->cols, EMPTY);
    
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            c->matrix[i][j] = a->matrix[i][j] - b->matrix[i][j];
        }
    }
    return c;
}
/* Multiply a matrix with a complex constant*/
Matrix *
consmul(pool         *p,
        Matrix       *a,
        complex double c)
        
{
    Matrix *b = matcreate(p, a->rows, a->cols, EMPTY);
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            b->matrix[i][j] = c * a->matrix[i][j];
        }
    }
    return b;
}
/*Copies column col1 from a to col2 in b*/
void 
matcopycol(Matrix       *a,
           int        col1,
           Matrix       *b,
           int        col2) 
{     
  for (int i = 0; i < a->rows; i++) {
    b->matrix[i][col2] = a->matrix[i][col1];
  }
}
/*Divide the matrix column "col" with the number "div"*/
Matrix *
matcoldiv(Matrix *a,
          int   col,
          double div)
{
    for(int i = 0; i < a->rows; i++) a->matrix[i][col]=a->matrix[i][col]/div;
    return a;
}
/*Subtracts matrix b's column col2 from matrix a's column col1*/
Matrix *
matcolsub(Matrix *a,
          int  col1,
          Matrix *b,
          int col2)
{
    for (int i = 0; i<a->cols; i++) {
        a->matrix[i][col1] -= b->matrix[i][col2];
    }
    return a;
}
           
/*Multiplies matrix column col with div*/
Matrix *
matcolmul(Matrix *a,
          int   col,
          double div)
          
{
    for(int i = 0; i < a->rows; i++) a->matrix[i][col]=a->matrix[i][col]*div;
    return a;
}
/*Produces a copy of matrix a and returns it.*/
Matrix *
matcopy (pool   *p,
         Matrix *a)
         
{
    Matrix *b = matcreate(p, a->rows, a->cols, EMPTY);
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            b->matrix[i][j] = a->matrix[i][j];
        }
    }
    return b;
}
/* print_matrix simply prints the matrix elements
 * row by row to standard output.
 */
void
print_matrix (Matrix *a)
{
    for(int i = 0; i < a->rows; i++) {
        for (int j = 0; j < a->cols; j++) {
            printf("%.2f %+.2fi  ", creal(a->matrix[i][j]),cimag(a->matrix[i][j]));
        }
        printf("\n");
    }
    printf("\n");
}

/* cmpmat compares to matrices by comparing each individual value.
 * if these values are approximately the same (tolerance 1e-4), then
 * it returns true. This function can't be used for very approximate comparisons
 * like when using infinite matrices.
 */
bool
cmpmat (Matrix *a,
        Matrix *b)
{
    if(a->cols != b->cols || a->rows != b->rows) {
        printf("Dimension error in matrix comparison");
        return false;
    } else {
        for(int i = 0; i < a->rows; i++) {
            for(int j = 0; j < a->cols; j++) {
                double ab = fabs(creal(a->matrix[i][j]) -\
                                         creal(b->matrix[i][j]));
                double ba = fabs(cimag(a->matrix[i][j]) -\
                                         cimag(b->matrix[i][j]));
                if(ab < 1e-4 && ba < 1e-4) continue;
                else return false;
            }
        }
        return true;
    }
}
/* Returns a matrix of the given type.
 * Annihilation matrix representation can be found from Wikipedia, for instance. Then creation
 * operator is just the dagger of that. That's why the code around caseANNIHILATION and case CREATIOn
 * is a bit awkward.
 */
Matrix *
matcreate (pool           *p,
           int            rows,
           int            cols,
           enum mat_types mat_type)
{   
    Matrix *a = pool_alloc(p, sizeof(Matrix));
    complex double **mat = pool_alloc(p, rows * sizeof(complex double *));
    switch(mat_type) {
    case EMPTY:
        for(int i = 0; i < rows; i++) {
            mat[i] = pool_alloc(p, cols * sizeof(complex double));
        }
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                mat[i][j] = 0.0;
            }
        }
        break;
    case ANNIHILATION: case CREATION: 
        for(int i = 0; i < rows; i++) {
            mat[i] = pool_alloc(p, cols * sizeof(complex double));
        }
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                double dj = (double) j;
                if (i == j-1) mat[i][j] = sqrt(dj);
                else mat[i][j] = 0.0;
            }
        }
        if(mat_type == CREATION) {
            a->matrix = mat;
            a->rows = rows;
            a->cols = cols;
            a = dagger(p, a);
        }
        break;

    case IDENTITY:
        for(int i = 0; i < rows; i++) {
            mat[i] = pool_alloc(p, cols * sizeof(complex double));
        }
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                if(i == j) mat[i][j] = 1.0;
                else mat[i][j] = 0.0;
            }
        }
        break;
    }
    if (mat_type != CREATION) {
        a->matrix = mat;
        a->cols = cols;
        a->rows = rows;
    }
    return a;
}
Matrix *
trim_zero(pool *p,
          Matrix *a)
{   
    int rowflag = 0;
    int colflag = 0;
    //int colflag = 0;
    int removedcols = 0;
    int removedrows = 0;
    //Checks if we find a 0-row at the beginning, break if not 
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            if(fabs(creal(a->matrix[i][j])) < 1e-3 && fabs(cimag(a->matrix[i][j])) < 1e-3) {
                ++rowflag;
            }
        }
        if(rowflag == a->cols) {
            //If we found a 0-row, check that the corresponding row on the other side of the matrix is also 0.
            //For example, if the first row is full of 0, the last one must also be, and if the second row has 0s,
            //so must the second to last, etc.
            rowflag = 0;
            
            for(int z = 0; z < a->cols; z++) {
                int k = a->rows - i - 1;
                if(fabs(creal(a->matrix[k][z])) < 1e-3 && fabs(cimag(a->matrix[k][z])) < 1e-3) {
                    ++rowflag;
                }
            }
            if(rowflag == a->cols) {
                ++removedrows;
                rowflag = 0;
            } else goto check;
            
        } else goto check;
    }
check:;    //Same logic as for above, this time checking for empty columns.
    for(int j = 0; j < a->cols; j++) {
        for(int i = 0; i < a->rows; i++) {
            if(fabs(creal(a->matrix[i][j])) < 1e-3 && fabs(cimag(a->matrix[i][j])) < 1e-3) {
                ++colflag;
            }
        }
        //If all numbers on the column are zero..
        if(colflag == a->rows) {
            colflag = 0;
            
            for(int z = 0; z < a->rows; z++) {
                int k = a->cols - j - 1;
                if(fabs(creal(a->matrix[z][k])) < 1e-3 && fabs(cimag(a->matrix[z][k])) < 1e-3) {
                    ++colflag;
                }
            }
            if(colflag == a->rows) {
                ++removedcols;
                colflag = 0;
            } else goto out;
        } else goto out;
    }
out:;
    Matrix *c = matcreate(p, a->rows - 2*removedrows, a->cols - 2*removedcols, EMPTY);
    for(int i = removedrows; i < a->rows-removedrows; i++) {
        for(int j = removedcols; j < a->cols-removedcols; j++) {
            c->matrix[i-removedrows][j-removedcols] = a->matrix[i][j];
        }
    }
    return c;
    
}
/* Basic matrix transpose routine */
Matrix *
dagger (pool *p, Matrix *a)
{

    Matrix *b = pool_alloc(p, sizeof(Matrix));
    complex double **mat = pool_alloc(p, a->cols * sizeof(complex double *));
    for(int i = 0; i < a->cols; i++) {
            mat[i] = pool_alloc(p, a->rows * sizeof(complex double));
    }
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
               mat[j][i] = conj(a->matrix[i][j]);
        }
    }
    b->rows = a->cols;
    b->cols = a->rows;
    b->matrix = mat;
    return b;
}
