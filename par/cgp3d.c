#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "omp.h"
#include "pcg.h"
#include "params.h"

/*@T
 * \section{3D Laplace operator}
 *
 * The 3D Laplacian looks like
 * \[
 *    (Ax)_{ijk} = h^{-2} 
 *    \left( 6x_{ijk} - \sum_{q,r,s: |q-i|+|j-r|+|k-s|=1} x_{qrs} \right).
 * \]
 * The [[mul_poisson3d]] function applies the 3D Laplacian to an
 * $n \times n \times n$ mesh of $[0,1]^3$ ($N = n^3$, $h = 1/(n-1)$),
 * assuming Dirichlet boundary conditions.
 *@c*/
void mul_poisson3d(int N, void* data, 
                   double* restrict Ax, 
                   double* restrict x)
{
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

    int n = *(int*) data;
    int inv_h2 = (n-1)*(n-1);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double xx = X(i,j,k);
                double xn = (i > 0)   ? X(i-1,j,k) : 0;
                double xs = (i < n-1) ? X(i+1,j,k) : 0;
                double xe = (j > 0)   ? X(i,j-1,k) : 0;
                double xw = (j < n-1) ? X(i,j+1,k) : 0;
                double xu = (k > 0)   ? X(i,j,k-1) : 0;
                double xd = (k < n-1) ? X(i,j,k+1) : 0;
                AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
            }
        }
    }

    #undef AX
    #undef X
}


/*@T
 * \section{Preconditioners for the Laplacian}
 *
 * \subsection{The identity preconditioner}
 *
 * The simplest possible preconditioner is the identity ([[pc_identity]]):
 *@c*/
void pc_identity(int n, void* data, double* Ax, double* x)
{
    memcpy(Ax, x, n*sizeof(double));
}

/*@T
 * \subsection{SSOR preconditioning}
 *
 * In terms of matrix splittings, if $A = L + D + L^T$ where $D$
 * is diagonal and $L$ is strictly lower triangular, the SSOR preconditioner
 * is given by
 * \[
 *   M(\omega) = \frac{1}{2-\omega} (D/\omega+L) (D/\omega)^{-1} (D/omega+L)^T,
 * \]
 * where $\omega$ is a relaxation parameter.  More algorithmically,
 * SSOR means looping through the unknowns and updating each by adding
 * $\omega$ times the Gauss-Seidel step; then doing the same thing,
 * but with the opposite order.  Choosing an optimal value of $\omega$
 * is not all that easy; there are heuristics when SOR is being used as
 * a stationary method, but I'm not sure that they apply when it is used
 * as a preconditioner.  The simplest thing to do is just to play with it.
 *
 * In order to apply SSOR preconditioning, we need the size $n$ of the
 * mesh (though in principle we could compute that from the number of
 * mesh points) as well as the parameter $\omega$.  We store these
 * parameters in a [[pc_ssor_p3d_t]] structure.
 *@c*/
typedef struct pc_ssor_p3d_t {
    int n;          /* Number of points in one direction on mesh */
    double omega;   /* SSOR relaxation parameter */
} pc_ssor_p3d_t;


/*@T
 *
 * The [[ssor_forward_sweep]], [[ssor_backward_sweep]], and [[ssor_diag_sweep]]
 * respectively apply $(D/\omega+L)^{-1}$, $(D/omega+L)^{-T}$, and
 * $(2-\omega) D/\omega$.  Note that we're okay with ignoring the $h^{-2}$
 * factor for computing the preconditioner --- multiplying $M$ by a scalar
 * constant doesn't change the preconditioned Krylov subspace.  Also note
 * that these functions can operate on just part of the mesh (rather than
 * the whole thing).  This will be useful shortly when we discuss additive
 * Schwarz preconditioners.
 *@c*/
void ssor_forward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                        double* restrict Ax, double w)
{
    #define min(a,b) ((a)<(b)?(a):(b))
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    #define BLOCK_SIZE 10        
    double xx, xn, xe, xu;
	
    // Inter block iteraions
    for (int d=0; d<= i2+j2+k2-3; d++) {
	for (int ii=i1;ii< min(d,i2-1);ii++)
	{
	int kk=k1;
	int jj = d-ii-kk;
	
		for (; jj>=j1 && kk<min(d,k2-1)-ii;jj--,kk++){
	//	printf("ii=%d jj=%d kk=%d d=%d\n",ii,jj,kk,d); 
		assert(d==ii+jj+kk);	
		{
		   for (int k = kk*BLOCK_SIZE; k < min(k2,kk*BLOCK_SIZE+BLOCK_SIZE); k++) {
                       for (int j = jj*BLOCK_SIZE; j < min(j2,jj*BLOCK_SIZE+BLOCK_SIZE); j++) {
                           for (int i = ii*BLOCK_SIZE; i < min(i2,ii*BLOCK_SIZE+BLOCK_SIZE); i++) {
                        	xx = AX(i,j,k);
                       		xn = (i > 0)   ? AX(i-1,j,k) : 0;
                       		xe = (j > 0)   ? AX(i,j-1,k) : 0;
                       		xu = (k > 0)   ? AX(i,j,k-1) : 0;
                       		AX(i,j,k) = (xx+xn+xe+xu)/6*w;
					}
 		                   }
		                }
           		 }

		}	
         }
 
	   
    #undef AX

}
}
void ssor_backward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                         double* restrict Ax, double w)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int k = k2-1; k >= k1; --k) {
        for (int j = j2-1; j >= j1; --j) {
            for (int i = i2-1; i >= i1; --i) {
                double xx = AX(i,j,k);
                double xs = (i < n-1) ? AX(i+1,j,k) : 0;
                double xw = (j < n-1) ? AX(i,j+1,k) : 0;
                double xd = (k < n-1) ? AX(i,j,k+1) : 0;
                AX(i,j,k) = (xx+xs+xw+xd)/6*w;
            }
        }
    }
    #undef AX
}

void ssor_diag_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                     double* restrict Ax, double w)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) *= (6*(2-w)/w);
    #undef AX
}

/* @T
 *
 * Finally, the [[pc_ssor_poisson3d]] function actually applies the
 * preconditioner.
 *@c*/
void pc_ssor_poisson3d(int N, void* data, 
                       double* restrict Ax, 
                       double* restrict x)
{
    pc_ssor_p3d_t* ssor_data = (pc_ssor_p3d_t*) data;
    int n = ssor_data->n;
    double w = ssor_data->omega;
    
    memcpy(Ax, x, N*sizeof(double));
    ssor_forward_sweep (n, 0,n, 0,n, 0,n, Ax, w);
    ssor_diag_sweep    (n, 0,n, 0,n, 0,n, Ax, w);
    ssor_backward_sweep(n, 0,n, 0,n, 0,n, Ax, w);
}

/*@T
 * \subsection{Additive Schwarz preconditioning}
 *
 * One way of thinking about Jacobi and Gauss-Seidel is as a sequence
 * of local relaxation operations, each of which updates a variable 
 * (or a set of variables, in the case of block variants) assuming
 * that the neighboring variables are known.  With Jacobi and Gauss-Seidel,
 * we update each variable exactly once in each pass.  In Schwarz methods,
 * we can update some variables with {\em multiple} relaxation operations
 * in a single pass.  To give an example, we will give an additive Schwarz
 * (Jacobi-like) method that updates the bottom half of the domain and
 * the top half of the domain with a little bit of overlap.  To update the
 * variables in each of the overlapping subdomains, we use one sweep of the
 * SSOR step described in the previous section.
 *
 * As mentioned in class, Schwarz-type preconditioners are a fantastic
 * match for distributed memory computation, since the processors only
 * communicate through the data in the overlap region.  Note, though, that the
 * specific variant I mentioned in class (restrictive additive Schwarz)
 * can't be used with conjugate gradient methods, because it does not yield
 * symmetric preconditioners.
 *
 * We describe the parameters for the Schwarz-type preconditioner with
 * SSOR-based approximate solves in a [[pc_schwarz_p3d_t]] structure.
 *@c*/
typedef struct pc_schwarz_p3d_t {
    int n;           /* Number of mesh points on a side */
    int overlap;     /* Number of points through the overlap region */
    double omega;    /* SSOR relaxation parameter */
    double* scratch; /* Scratch space used by the preconditioner */
} pc_schwarz_p3d_t;


/*@T
 *
 * In order to compute independently on overlapping subdomains, we first
 * get the local data from the vector to which we're applying the
 * preconditioner; then we do an inexact solve on the local piece of
 * the data; and then we write back the updates from the solve.
 * The data motion is implemented in [[schwarz_get]] and [[schwarz_add]].
 *@c*/
void schwarz_get(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                 double* restrict x_local, 
                 double* restrict x)
{
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define XL(i,j,k) (x_local[((k)*n+(j))*n+(i)])

    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k) = X(i,j,k);

    if (k1 > 0)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k1-1) = 0;
    if (j1 > 0)
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j1-1,k) = 0;
    if (i1 > 0)
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i1-1,j,k) = 0;

    if (k2 < n)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k2) = 0;
    if (j2 < n)
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j2,k) = 0;
    if (i2 < n)
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i2,j,k) = 0;

    #undef XL
    #undef X
}

void schwarz_add(int n, int i1, int i2, int j1, int j2, int k1, int k2, 
                 double* restrict Ax_local, 
                 double* restrict Ax)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    #define AXL(i,j,k) (Ax_local[((k)*n+(j))*n+(i)])
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) += AXL(i,j,k);
    #undef AXL
    #undef AX
}

/*@T
 *
 * The [[pc_schwarz_poisson3d]] function applies a preconditioner by
 * combining independent SSOR updates for the bottom half plus an overlap
 * region (a slap $n/2+o/2$ nodes thick), then updating the top half plus
 * an overlap region (another $n/2+o/2$ node slab).  The same idea could
 * be applied to more regions, or to better approximate solvers.
 *@c*/
void pc_schwarz_poisson3d(int N, void* data, 
                          double* restrict Ax, 
                          double* restrict x)
{
    pc_schwarz_p3d_t* ssor_data = (pc_schwarz_p3d_t*) data;
    double* scratch = ssor_data->scratch;
    int n = ssor_data->n;
    int o = ssor_data->overlap/2;
    double w = ssor_data->omega;
    memset(Ax, 0, N*sizeof(double));

    int n1 = n/2+o;
    int n2 = n/2-o;

    schwarz_get        (n, 0,n, 0,n, 0,n1, scratch, x);
    ssor_forward_sweep (n, 0,n, 0,n, 0,n1, scratch, w);
    ssor_diag_sweep    (n, 0,n, 0,n, 0,n1, scratch, w);
    ssor_backward_sweep(n, 0,n, 0,n, 0,n1, scratch, w);
    schwarz_add        (n, 0,n, 0,n, 0,n1, scratch,Ax);

    schwarz_get        (n, 0,n, 0,n, n2,n, scratch, x);
    ssor_forward_sweep (n, 0,n, 0,n, n2,n, scratch, w);
    ssor_diag_sweep    (n, 0,n, 0,n, n2,n, scratch, w);
    ssor_backward_sweep(n, 0,n, 0,n, n2,n, scratch, w);
    schwarz_add        (n, 0,n, 0,n, n2,n, scratch,Ax);
}

/*@T
 * \section{Forcing functions}
 *
 * The convergence of CG depends not only on the operator and the
 * preconditioner, but also on the right hand side.  If the error
 * is very high frequency, the convergence will appear relatively
 * fast.  Without a good preconditioner, it takes more iterations
 * to correct a smooth error.  In order to illustrate these behaviors,
 * we provide two right-hand sides: a vector with one nonzero 
 * (computed via [[setup_rhs0]]) and a vector corresponding to a smooth
 * product of quadratics in each coordinate direction ([[setup_rhs1]]).
 *@c*/
void setup_rhs0(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    b[0] = 1;
}

void setup_rhs1(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    for (int i = 0; i < n; ++i) {
        double x = 1.0*(i+1)/(n+1);
        for (int j = 0; j < n; ++j) {
            double y = 1.0*(i+1)/(n+1);
            for (int k = 0; k < n; ++k) {
                double z = 1.0*(i+1)/(n+1);
                b[(k*n+j)*n+i] = x*(1-x) * y*(1-y) * z*(1-z);
            }
        }
    }
}

/*@T
 * \section{The [[main]] event}
 *
 * The main driver is pretty simple: read the problem and solver parameters
 * using [[get_params]] and then run the preconditioned solve.
 *
 * I originally thought I'd write your own option routine, but I changed my
 * mind!
 *@c*/
int main(int argc, char** argv)
{
    solve_param_t params;
    if (get_params(argc, argv, &params))
        return -1;

    int n = params.n;
    int N = n*n*n;

    double* b = malloc(N*sizeof(double));
    double* x = malloc(N*sizeof(double));
    double* r = malloc(N*sizeof(double));
    memset(b, 0, N*sizeof(double));
    memset(x, 0, N*sizeof(double));
    memset(r, 0, N*sizeof(double));

    /* Set up right hand side */
    setup_rhs1(n,b);

    /* Solve via PCG */
    int maxit   = params.maxit;
    double rtol = params.rtol;

    if (params.ptype == PC_SCHWARZ) {
        double* scratch = malloc(N*sizeof(double));
        pc_schwarz_p3d_t pcdata = {n, params.overlap, params.omega, scratch};
        pcg(N, pc_schwarz_poisson3d, &pcdata, mul_poisson3d, &n, x, b, 
            maxit, rtol);
        free(scratch);
    } else if (params.ptype == PC_SSOR) {
        pc_ssor_p3d_t ssor_data = {n, params.omega};
        pcg(N, pc_ssor_poisson3d, &ssor_data, mul_poisson3d, &n, x, b, 
            maxit, rtol);
    } else {
        pcg(N, pc_identity, NULL, mul_poisson3d, &n, x, b, maxit, rtol);
    }

    /* Check answer */
    mul_poisson3d(N, &n, r, x);
    double rnorm2 = 0;
    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];
    for (int i = 0; i < n; ++i) rnorm2 += r[i]*r[i];
    printf("rnorm = %g\n", sqrt(rnorm2));

    free(r);
    free(x);
    free(b);
}
