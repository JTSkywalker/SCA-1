#ifndef PARAMS_H
#define PARAMS_H

/*@T
 * \section{Solver parameters}
 * 
 * The [[solve_param_t]] structure holds the parameters that
 * describe the simulation.  These parameters are filled in
 * by the [[get_params]] function.  Details of the parameters
 * are described elsewhere in the code.
 *@c*/
enum {                /* Types of preconditioners available: */
    PC_ID = 1,        /* 1. Identity                         */
    PC_SSOR = 2,      /* 2. SSOR                             */
    PC_SCHWARZ = 3    /* 3. Additive Schwarz                 */
};
    
typedef struct solve_param_t {
    int    n;       /* Mesh size */
    int    maxit;   /* Maximum PCG iterations */
    double rtol;    /* Relative residual convergence tolerance */
    int    ptype;   /* Preconditioner type */
    double omega;   /* SSOR relaxation parameter */
    int    overlap; /* Overlap size */
} solve_param_t;

int get_params(int argc, char** argv, solve_param_t* params);

/*@q*/
#endif /* PARAMS_H */
