#ifndef PCG_H
#define PCG_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef void (*mul_fun_t)(int n, void* data, double* Ax, double* x);
double dot(int n, const double* x, const double* y);
double pcg(int n,
           mul_fun_t Mfun, void* Mdata,
           mul_fun_t Afun, void* Adata,
           double* restrict x,
           const double* restrict b,
           int maxit,
           double rtol);

#endif /* PCG_H */
