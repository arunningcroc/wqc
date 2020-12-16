#include <stdio.h>
#include <stdlib.h>
#include "operators.h"
#include "pool.h"
#include "ode.h"
#include <math.h>
/* This file is only for testing the code. You can run the program
 * by opening the exe named scquantumoptics, and then
 * that will run these tests. But for actually using the code, one should
 * use the Python bindings.
 */
int
test_operators ()
{
    int sum = 0;
    enum mat_types t = EMPTY;
    pool *p = begin_calculation(10000*sizeof(double));
    //Test matrix multiplication
    Matrix *a = matcreate(p, 2, 2, t);
    Matrix *b = matcreate(p, 2, 2, t);
    //Test matrix 1
    a->matrix[0][0] = 2;
    a->matrix[0][1] = 3+2I;
    a->matrix[1][0] = 4;
    a->matrix[1][1] = 5;
    //Test matrix 2
    b->matrix[0][0] = 5+3I;
    b->matrix[0][1] = 1;
    b->matrix[1][0] = 6;
    b->matrix[1][1] = 7;
    Matrix *c = matmul(p, a,b);
    Matrix *d = matcreate(p, 2, 2, t);
    d->matrix[0][0] = 28+18I;
    d->matrix[0][1] = 23+14I;
    d->matrix[1][0] = 50+12I;
    d->matrix[1][1] = 39+0I;
    if(cmpmat(c,d)) sum += 1;
    //Test transpose for block matrix
    Matrix *e = dagger(p,d);
    Matrix *dtransp = matcreate(p,2,2,t);

    dtransp->matrix[0][0] = 28-18I;
    dtransp->matrix[0][1] = 50-12I;
    dtransp->matrix[1][0] = 23-14I;
    dtransp->matrix[1][1] = 39-0.0I;
    if(cmpmat(e,dtransp)) sum += 1;

    //Test transposing a non-block matrix
    Matrix *f = matcreate(p, 1, 2, t);
    f->matrix[0][0] = 1;
    f->matrix[0][1] = 2;
    Matrix *ftransp = matcreate(p, 2, 1, t);
    ftransp->matrix[0][0] = 1;
    ftransp->matrix[1][0] = 2;
    Matrix *df = dagger(p,f);
    if(cmpmat(df,ftransp)) sum += 1;
    //Tets the Kronecker delta (tensor)
    Matrix *tenstest = matcreate(p, 2, 2, EMPTY);
    Matrix *tenstest2 = matcreate(p, 2, 2, EMPTY);
    tenstest->matrix[0][0] = 1;
    tenstest->matrix[0][1] = 2;
    tenstest->matrix[1][0] = 3;
    tenstest->matrix[1][1] = 4;

    tenstest2->matrix[0][0] = 0;
    tenstest2->matrix[0][1] = 5;
    tenstest2->matrix[1][0] = 6;
    tenstest2->matrix[1][1] = 7;

    Matrix *tensres = tensor(p, tenstest, tenstest2);
    Matrix *realres = matcreate(p, 4, 4, EMPTY);
    realres->matrix[0][0] = 0;
    realres->matrix[0][1] = 5;
    realres->matrix[0][2] = 0;
    realres->matrix[0][3] = 10;
    
    realres->matrix[1][0] = 6;
    realres->matrix[1][1] = 7;
    realres->matrix[1][2] = 12;
    realres->matrix[1][3] = 14;
    
    realres->matrix[2][0] = 0;
    realres->matrix[2][1] = 15;
    realres->matrix[2][2] = 0;
    realres->matrix[2][3] = 20;
    
    realres->matrix[3][0] = 18;
    realres->matrix[3][1] = 21;
    realres->matrix[3][2] = 24;
    realres->matrix[3][3] = 28;

    if(cmpmat(tensres,realres)) sum += 1;
    Matrix *identity = matcreate(p, 2, 2, IDENTITY);
    Matrix *idreal = matcreate(p, 2, 2, EMPTY);
    idreal->matrix[0][0] = 1.0;
    idreal->matrix[0][1] = 0.0;
    idreal->matrix[1][0] = 0.0;
    idreal->matrix[1][1] = 1.0;
    if(cmpmat(identity,idreal)) sum += 1;
    if(trace(identity) == 2) sum += 1;
    
    Matrix *testeig = matcreate(p, 4, 4, EMPTY);
   
    
    testeig->matrix[0][0] = 0;
    testeig->matrix[0][1] = 0;
    testeig->matrix[0][2] = 0;
    testeig->matrix[0][3] = 0;
    
    testeig->matrix[1][0] = 0;
    testeig->matrix[1][1] = 2;
    testeig->matrix[1][2] = 3;
    testeig->matrix[1][3] = 0;
    
    testeig->matrix[2][0] = 0;
    testeig->matrix[2][1] = 4;
    testeig->matrix[2][2] = 5;
    testeig->matrix[2][3] = 0;
    
    testeig->matrix[3][0] = 0;
    testeig->matrix[3][1] = 0;
    testeig->matrix[3][2] = 0;
    testeig->matrix[3][3] = 0;
    

    //print_matrix(testeig);
    //complex double *eigenvalues = solve_eigenvalues(p, testeig);
    //Matrix *trimmed = trim_zero(p, testeig);
    //print_matrix(trimmed);
    printf("\n");
    //int number = 25;
    /*Test more eigenvalues*/
    /*Matrix *annt = matcreate(p, number, number, ANNIHILATION);
    Matrix *ident = matcreate(p, number, number, IDENTITY);
    
    Matrix *a1 = tensor(p, annt, ident);
    //print_matrix(a1);
    Matrix *a2 = tensor(p, ident, annt);
    Matrix *ad1 = dagger(p, a1);
    Matrix *ad2 = dagger(p, a2);
    Matrix *numberop = matsum(p, matmul(p, ad1, a1), matmul(p, ad2, a2));
    Matrix *shift = consmul(p,tensor(p,ident, ident), -number);
    Matrix *totalshift = matsum(p, numberop, shift);
    Matrix *Heps = (matsub(p, matmul(p, ad1, a1), matmul(p, ad2, a2)));
    Matrix *Hv   = (matsum(p, matmul(p, ad1, a2), matmul(p, ad2, a1)));
    Matrix *Hc   = consmul(p,(matsub(p, matmul(p, ad1, a1), matmul(p, ad2, a2))),0.5);
    Matrix *Hcsqrd = matmul(p, Hc, Hc);
    Matrix *Htotal = matsum(p, Heps,matsum(p, Hv, Hcsqrd));
    Matrix *Hshiftterm = matsub(p, Htotal, consmul(p, totalshift, 10000));
    Matrix *Htrim = trim_zero(p, Hshiftterm);
    Matrix **dimer = solve_eigenvalues(p, Htrim);
    complex double *eigvals = pool_alloc(p, dimer[1]->rows * sizeof(complex double));
    for(int i = 0; i < dimer[1]->rows; i++) {
        if(fabs(creal(dimer[1]->matrix[i][i])) < 10*number)
            eigvals[i] = dimer[1]->matrix[i][i];
    }
    double biggest = -500;
    for(int i = 0; i < dimer[1]->rows; i++) {
        if(creal(eigvals[i]) > biggest) {
            biggest = eigvals[i];
        }
    }
    printf("%f", biggest);*/
    end_calculation(p);
    pool *new = begin_calculation(1000000*sizeof(double));
    double *jceig = pool_alloc(new, 200 * sizeof(double));
    double n = 0.0;
    double omega = 1;
    double omega0 = 1;
    double lambda = 1;
    Matrix *HJC = matcreate(new, 2, 2, EMPTY);
    Matrix **dimss;
    double bigst = -500;
    for(int i = 0; i < 100; i++) {
        HJC->matrix[0][0] = i*omega;
        HJC->matrix[0][1] = lambda * sqrt(i+1);
        HJC->matrix[1][0] = lambda * sqrt(i+1);
        HJC->matrix[1][1] = (i+1)*omega - 0.5*omega0;
        dimss = solve_eigenvalues(new, HJC);
        for(int j = 0; j < 2; j++) {
            jceig[i+j] = dimss[1]->matrix[j][j];
        }
        
    }
    FILE *nss = fopen("eigsjc.txt", "w");
    for(int i = 0; i < 100; i++) {
        fprintf(nss, "%d %f\n", i, jceig[i]);
        fprintf(nss, "%d %f\n", i, jceig[i+1]);
    }
    fclose(nss);
    end_calculation(new);
    solve_bhlind();
    if(sum == 6) return 1;
    else return 0;

}
int test_ode() {
    enum prob_types o = TEST;
    solve_ode(o, 1,0.5, 0.5, 2 ,1);

    return 1;
}
int
main (int argc,
      char *argv[])
{
    if (test_operators()) printf("Operator functions OK\n");
    else printf("Operator function tests failed. Rebuild\n");
    if (test_ode()) printf("ODE solvers OK\n");
    else printf("ODE solver tests failed. Rebuild.\n");

    return EXIT_SUCCESS;
}
