/* Routine for solving ODEs. Here it is imlemented as just an ordinary RK4. */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "ode.h"
#include "operators.h"

complex double
testf(double x, double t)
{
    return -pow(x,3) + sin(t);
}
complex double *
rabi(pool           *p,
     complex double  x,
     complex double x2,
     double          t,
     Matrix    **basis,
     Matrix         *d,
     double  fstrength,
     double    detuning)
{
    //d is the dipole moment matrix, which is offdiagonal and depends mostly on field strength.
    d->matrix[0][1] = fstrength;
    d->matrix[1][0] = fstrength;
    Matrix *b0 = dagger(p, basis[0]);
    Matrix *b1 = dagger(p, basis[1]);
    Matrix *k0 = basis[0];
    Matrix *k1 = basis[1];
    //Eq 4.20 from the book.
    complex double C1 = I/2 *matmul(p,b0,matmul(p,d,k1))->matrix[0][0]*x2*cexp(I*detuning*t);
    complex double C2 = I/2 *matmul(p,b1,matmul(p,d,k0))->matrix[0][0]*x*cexp(-I*detuning*t);
    complex double *sols = pool_alloc(p, 2*sizeof(complex double));
    sols[0] = C1;
    sols[1] = C2;
    return sols;
    
}

complex double *
rabi_norwa(pool           *p,
           complex double  x,
           complex double x2,
           double          t,
           Matrix    **basis,
           Matrix         *d,
           double     wfield,
           double      watom,
           double fstrength)
{
    //d is the dipole moment matrix, which is offdiagonal and depends mostly on field strength.
    d->matrix[0][1] = fstrength;
    d->matrix[1][0] = fstrength;
    Matrix *b0 = dagger(p, basis[0]);
    Matrix *b1 = dagger(p, basis[1]);
    Matrix *k0 = basis[0];
    Matrix *k1 = basis[1];
    //Eq 4.20 from the book
    complex double C1 = I/2 *matmul(p,b0,matmul(p,d,k1))->matrix[0][0]*x2*cos(wfield*t)*cexp(I*watom*t);
    complex double C2 = I/2 *matmul(p,b1,matmul(p,d,k0))->matrix[0][0]*x*cos(wfield*t)*cexp(-I*watom*t);
    complex double *sols = pool_alloc(p, 2*sizeof(complex double));
    sols[0] = C1;
    sols[1] = C2;
    return sols;
    
}
void
solve_ode(enum prob_types typ,
          double       wfield,
          double        watom,
          double    fstrength,
          int               n,
          int             iter)
{
    //Allocate excessive memory because nobody cares about kilobytes.
    pool *p = begin_calculation(2000000*sizeof(double));
    double detuning = wfield - watom;
    //Integration parameters (can't change for now - gets the point across)
    int niter = 10000;
    double a = 0.0;
    double b = 25.0;
    double N = (double)niter;
    double h = (b-a)/N;
    
    //-------
    char name[12];
    sprintf(name, "data%d.txt", iter);
    
    //Used for holding the actual values of the integration, so basically
    //this is C(t). There are separate arrays for each individual equation,
    //so for each C(t).
    complex double **values = pool_alloc(p, n*sizeof(complex double *));
    for(int i = 0; i < n; i++) {
        values[i] = pool_alloc(p, niter*sizeof(complex double)); 
    }

    double *t = pool_alloc(p, niter*sizeof(double));
    complex double *xs = pool_alloc(p, niter*sizeof(complex double));
    for(int i = 0; i < n; i++) {
        if (i == 0) xs[i] = 1;
        else xs[i] = 0;
    }
    Matrix **altbasis = fock_basis(p, n);
    //E-field
    Matrix *d = matcreate(p, n, n, EMPTY);
    switch(typ){
    case TEST:
        for(int i = 0; i < 100; i++) {
            t[i] = 0.0 + i*0.1;
        }
        for(int i = 0; i < 100; i++) {
            values[0][i] = xs[0];
            double k1 = h*testf(xs[0],t[i]);
            double k2 = h*testf(xs[0]+0.5*k1, t[i]+0.5*h);
            xs[0] += k2;
        }
        FILE *f = fopen("result.txt", "w");
        if (f == NULL) printf("File error");

        for(int i=0; i < 100; i++) {
            fprintf(f, "%f %f\n", t[i], creal(values[0][i]));
        }
        fclose(f);
        end_calculation(p);
        break;
    //Case for the 2-level atom Rabi model under RWA approximation
    case RABI:
        for(int i = 0; i < niter; i++) {
            t[i] = 0.0 + h*i;
        }
        for(int i = 0; i < niter; i++) {
            //This is standard euler integration
            values[0][i] = xs[0];
            values[1][i] = xs[1];
            complex double *sols = rabi(p, xs[0], xs[1], t[i], altbasis, d, fstrength,detuning);
            xs[0] += h*sols[0];
            xs[1] += h*sols[1]; 

        }
        FILE *f2 = fopen(name, "w");
        if (f2 == NULL) printf("File error");

        for(int i=0; i < niter; i++) {
            double cc = conj(values[1][i]) * values[1][i];
            fprintf(f2, "%f %f\n", t[i], cc);
        }
        end_calculation(p);
        break;
    //Rabi with no RWA
    case NORWA:
        for(int i = 0; i < niter; i++) {
            t[i] = 0.0 + h*i;
        }
        for(int i = 0; i < niter; i++) {
        //Euler integration
            values[0][i] = xs[0];
            values[1][i] = xs[1];
            complex double *sols = rabi_norwa(p, xs[0], xs[1], t[i], altbasis, d, wfield, watom, fstrength);
            xs[0] += h*sols[0];
            xs[1] += h*sols[1]; 

        }
        FILE *f3 = fopen(name, "w");
        if (f3 == NULL) printf("File error");

        for(int i=0; i < niter; i++) {
            double cc = conj(values[1][i]) * values[1][i];
            fprintf(f3, "%f %f\n", t[i], cc);
        }
        end_calculation(p);
        break;
    }
        
}
Matrix *
bh_lindblad_form(pool             *p,
                 Matrix          *a1,
                 Matrix         *ad1,
                 Matrix          *a2,
                 Matrix         *ad2,
                 Matrix           *H,
                 Matrix         *rho,
                 complex double gamma)
{
    Matrix *commutator = matsub(p, matmul(p, H, rho), matmul(p, rho, H));
    Matrix *icomm = consmul(p, commutator, -I);
    
    Matrix *part1 = matmul(p, ad2, matmul(p, a2, rho));
    Matrix *part2 = matmul(p, rho, matmul(p, ad2, a2));
    Matrix *part3_noconst = matmul(p, a2, matmul(p, rho, ad2));
    Matrix *part3 = consmul(p, part3_noconst, -2.0);
    Matrix *brackets = matsum(p, part1, matsum(p, part2, part3));
    Matrix *brackets_with_multiplier = consmul(p, brackets, -gamma/2);
    
    Matrix *newrho = matsum(p, icomm, brackets_with_multiplier);
    return newrho;
}
void
solve_bhlind()
{
    int N = 2;
    int Np = N+1;
    int niter = 4000;
    pool *p = begin_calculation(300000000 * sizeof(double));
    
    Matrix *annt = matcreate(p, Np, Np, ANNIHILATION);
    Matrix *ident = matcreate(p, Np, Np, IDENTITY);
    
    Matrix *a1 = tensor(p, annt, ident);
    //print_matrix(a1);
    //Elements of the Hamiltonian
    Matrix *a2 = tensor(p, ident, annt);
    Matrix *ad1 = dagger(p, a1);
    Matrix *ad2 = dagger(p, a2);
    Matrix *Hv   = consmul(p, (matsum(p, matmul(p, ad1, a2), matmul(p, ad2, a1))), 0.3);
    Matrix *Hc   = consmul(p,(matsub(p, matmul(p, ad1, a1), matmul(p, ad2, a2))),(0.6/0.5));
    Matrix *Hcsqrd = matmul(p, Hc, Hc);
    Matrix *Htotal = matsum(p, Hv, Hcsqrd);
    
    complex double gamma = 0.02;
    Matrix *psi0_1 = matcreate(p, Np, 1, EMPTY);
    psi0_1->matrix[0][0] = 1;
    
    Matrix *psi0_2 = matcreate(p, Np, 1, EMPTY);
    psi0_2->matrix[Np-1][0] = 1;
    //Total wave function PSI of the state
    Matrix *psi0 = tensor(p, psi0_1, psi0_2);
    //Now, the initial density matrix.
    Matrix *rho = matmul(p, psi0, dagger(p,psi0));
    print_matrix(rho);
    double nk = 10.0;
    double dt = 0.05;
    double *t = pool_alloc(p, niter*sizeof(double));
    for(int i = 0; i < niter; i++) {
        t[i] = 0.0 + i * dt;
    }
    double *n1 = pool_alloc(p, niter*sizeof(double));
    for(int i = 0; i < niter; i++) {
        n1[i] = 0.0;
    }
    double *n2 = pool_alloc(p, niter*sizeof(double));
    for(int i = 0; i < niter; i++) {
        n2[i] = 0.0;
    }
    for(int i = 0; i < niter; i++) {
        printf("%d\n", i);
        n2[i] = creal(trace(matmul(p,rho, matmul(p, ad2, a2))));
        n1[i] = creal(trace(matmul(p,rho, matmul(p, ad1, a1))));
        for(int k = 0; k < nk; k++) {
            Matrix *new = matsum(p, rho, consmul(p,bh_lindblad_form(p, a1, ad1, a2, ad2, Htotal, rho, gamma),dt/nk));
            Matrix *dif = consmul(p, matsum(p, rho, new), 0.5);
            rho = matsum(p,rho, consmul(p,bh_lindblad_form(p, a1, ad1, a2, ad2, Htotal, dif, gamma),dt/nk));
        }
    }
    FILE *f3 = fopen("testi.txt", "w");
    for(int i=0; i < niter; i++) {
            fprintf(f3, "%f %f\n", t[i], n2[i]/N);
        }
    FILE *f4 = fopen("tiheys.txt", "w");
    for(int i=0; i < Np; i++) {
		fprintf(f4,"%f %f\n", t[i], creal(rho->matrix[i][i]));
		

    }
    fclose(f4);
    fclose(f3);
    end_calculation(p);
}
