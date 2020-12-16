#ifndef ODE
#define ODE
#include "pool.h"
enum prob_types {
    TEST, RABI, NORWA
};
void solve_ode(enum prob_types typ, double wfield, double watom, double fstrength, int n, int iter);
void solve_bhlind();
#endif
