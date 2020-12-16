from libcpp cimport bool
from libcpp cimport complex
cdef extern from "pool.h":
    ctypedef struct pool:
        char *next
        char *end
    pool *begin_calculation(size_t size)
    void end_calculation(pool *p)
    size_t pool_available(pool *p)
    void *pool_alloc(pool *p, size_t size)
cdef extern from "operators.h":
    ctypedef struct Matrix:
        double complex **matrix
        int rows
        int cols
    cdef enum mat_types:
        EMPTY = 0
        CREATION = 1
        ANNIHILATION = 2
    Matrix *matmul (pool *p, Matrix *a, Matrix *b)
    void print_matrix(Matrix *a)
    bool cmpmat(Matrix *a, Matrix *b)
    Matrix *matcreate(pool *p, int rows, int cols, mat_types mat_type)
    Matrix *dagger (pool *p, Matrix *a)
    Matrix **fock_basis (pool *p, int n)
    double vnorm(Matrix *a, int col)
    Matrix *matcoldiv(Matrix *a, int col, double div)
    Matrix *matcolsub(Matrix *a, int col1, Matrix *b, int col2)
    Matrix *matcolmul(Matrix *a, int col, double div)
    Matrix *matcopy (pool *p, Matrix *a)
cdef extern from "ode.h":
    cdef enum prob_types:
        TEST = 0
        RABI = 1
        NORWA = 2
        IDENTITY = 3
    void solve_ode(prob_types typ, double wfield, double watom, double fstrength, int n, int iter)
