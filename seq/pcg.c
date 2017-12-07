#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "timing.h"
#include "pcg.h"

#include <scorep/SCOREP_User.h>


double dot(int n, const double* x, const double* y)
{
    double result = 0;
    for (int i = 0; i < n; ++i)
        result += x[i]*y[i];
    return result;
}

/*@T
 * \section{Preconditioned CG}
 *
 * The PCG routine multiplies by $A$ and $M^{-1}$ through the
 * [[Mfun]] and [[Afun]] function pointers (taking a vector length $n$,
 * an opaque data object, an output buffer, and an input vector as arguments).
 * We also pass [[Mdata]] and [[Adata]] as a way of getting context into
 * the function calls%
 * \footnote{This could admittedly be more convenient in C}.
 * In addition, we take storage for the solution (set to an initial guess
 * on input) and the right hand side, as well as the maximum number of
 * iterations allowed and a relative error tolerance.
 *
 * The relative error tolerance is actually slightly subtle; we terminate
 * the iteration when
 * \[
 *   \frac{\|r^{(k)}\|_{M^{-1}}}
 *        {\|r^{(0)}\|_{M^{-1}}} < \mathrm{tol},
 * \]
 * where $\|\cdot\|_{M^{-1}}$ refers to the norm induced by the $M^{-1}$
 * inner product, i.e. $\|z\|_{M^{-1}}^2 = z^T M^{-1} z$.  This may or
 * may not be the norm anyone actually cares about... but it surely is
 * cheap to compute.
 *@c*/
double pcg(int n,
           mul_fun_t Mfun, void* Mdata,
           mul_fun_t Afun, void* Adata,
           double* restrict x,
           const double* restrict b,
           int maxit,
           double rtol)
{
    double* r = malloc(n*sizeof(double));
    double* z = malloc(n*sizeof(double));
    double* q = malloc(n*sizeof(double));
    double* p = malloc(n*sizeof(double));

    double rho0     = 0;
    double rho      = 0;
    double rho_prev = 0;
    double rtol2 = rtol*rtol;
    int is_converged = 0;
    int step;

    tic(0);

    /* Form residual */
    Afun(n, Adata, r, x);
    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];

    for (step = 0; step < maxit && !is_converged; ++step) {
        Mfun(n, Mdata, z, r);
        rho_prev = rho;
        rho = dot(n, r, z);
        if (step == 0) {
            rho0 = rho;
            memcpy(p, z, n*sizeof(double));
        } else {
            double beta = rho/rho_prev;
	    SCOREP_USER_REGION_DEFINE(forP)
	      SCOREP_USER_REGION_BEGIN(forP, "forP", SCOREP_USER_REGION_TYPE_COMMON)
            for (int i = 0; i < n; ++i) p[i] = z[i] + beta*p[i];
	    SCOREP_USER_REGION_END(forP)
        }
        Afun(n, Adata, q, p);
        double alpha = rho/dot(n, p, q);
	SCOREP_USER_REGION_DEFINE(forPQ)
	  SCOREP_USER_REGION_BEGIN(forPQ, "forPQ", SCOREP_USER_REGION_TYPE_COMMON)
	  for (int i = 0; i < n; ++i) {
	    x[i] += alpha*p[i];
	    r[i] -= alpha*q[i];
	  }
	SCOREP_USER_REGION_END(forPQ)
        is_converged = (rho/rho0 < rtol2);
    }

    printf("%d steps, residual reduction %g (%s tol %g); time %g\n",
           step, sqrt(rho/rho0), is_converged ? "<=" : ">", rtol, toc(0));

    free(p);
    free(q);
    free(z);
    free(r);

    return rho/rho0;
}
