#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "timing.h"
#include "pcg.h"

//#include <scorep/SCOREP_User.h>


double dot(int n, const double* x, const double* y)
{
    double result = 0;    
#pragma omp parallel for schedule(guided) reduction(+:result) firstprivate(x,y)
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
    double t=0, time_dot=0;

    tic(0);
    tic(1);
    /* Form residual */
    Afun(n, Adata, r, x);
t+=toc(1);
    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];

    for (step = 0; step < maxit && !is_converged; ++step) {
      
        Mfun(n, Mdata, z, r);
        rho_prev = rho;
        tic(2);
	rho = dot(n, r, z);
        time_dot+=toc(2);	
	if (step == 0) {
            rho0 = rho;
            memcpy(p, z, n*sizeof(double));
        } else {
            double beta = rho/rho_prev;
	    // scorep instrumentation code
//	    SCOREP_USER_REGION_DEFINE(forP)
//	      SCOREP_USER_REGION_BEGIN(forP, "forP", SCOREP_USER_REGION_TYPE_COMMON)

	      #pragma omp parallel for
	      for (int i=0; i< n; i++)
		p[i] = z[i] + beta*p[i];
	   
//	    SCOREP_USER_REGION_END(forP)
        }
	tic(1);
        Afun(n, Adata, q, p);
	t+=toc(1);
	tic(2);
        double alpha =dot(n, p, q);
time_dot+=toc(2);
alpha=rho/alpha;
//	SCOREP_USER_REGION_DEFINE(forPQ)
//	  SCOREP_USER_REGION_BEGIN(forPQ, "forPQ", SCOREP_USER_REGION_TYPE_COMMON)
	  // It doesn't change ANYTHING! WHY ?!
	  //	#pragma omp parallel for schedule(static) firstprivate(p,q)
	  #pragma omp parallel for
	  for (int i = 0; i < n; ++i) {
	    x[i] += alpha*p[i];
	    r[i] -= alpha*q[i];
	  }
	
//	SCOREP_USER_REGION_END(forPQ)
        is_converged = (rho/rho0 < rtol2);
    }

    printf("%d steps, residual reduction %g (%s tol %g); time %g\n",
           step, sqrt(rho/rho0), is_converged ? "<=" : ">", rtol, toc(0));

    free(p);
    free(q);
    free(z);
    free(r);
    printf("time doing Afun=%lf\n", t);
    printf("time doing dot products=%lf\n", time_dot);
    return rho/rho0;
}
