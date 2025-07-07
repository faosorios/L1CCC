/* ID: RNG.h, last updated 2025-06-26, F.Osorio */

#ifndef RNG_H
#define RNG_H

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

/* some definitions */
#ifndef FCONE
# define FCONE
#endif
#define CUBE(x)         R_pow_di(x, 3)
#define DNULLP          (double *) 0
#define EPS_CONV        1.0e-2
#define FOURTH(x)       R_pow_di(x, 4)
#define GOLDEN          0.3819660112501051
#define IZERO(x)        (((x) == 0) ? 1 : 0)
#define MAX(a,b)        (((a)>(b)) ? (a) : (b))
#define MIN(a,b)        (((a)<(b)) ? (a) : (b))
#define OFFSET(n, inc)  (((inc) > 0) ? 0 : ((n) - 1) * (-(inc)))
#define repeat          for(;;)
#define SGN(x)          (((x) >= 0) ? 1.0 : -1.0)
#define SQR(x)          R_pow_di(x, 2)

/* multivariate random number generators */
void RNG_cauchy(double *, int *, int *, double *, double *);
void RNG_contaminated(double *, int *, int *, double *, double *, double *, double *);

/* matrix computations */
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void chol_decomp(double *, int, int, int, int *);
extern void mult_triangular_mat(double, double *, int, int, int, char *, char *, char *, char *, double *, int);
extern void scale_inplace(double *, int, int, double);

/*DEBUG */
extern void print_mat(double *, int, int, int, char *);

#endif /* RNG_H */
