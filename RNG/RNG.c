/* ID: RNG.c, last updated 2025-06-26, F.Osorio */

#include "RNG.h"

/* static functions.. */
static void rmcauchy_std(double *, int, int);
static void rmcontaminated_std(double *, double, double, int, int);
/* ..end declarations */

/* ========================================================================== *
 * multivariate Cauchy random generation
 * ========================================================================== */

void
RNG_cauchy(double *y, int *nobs, int *nvar, double *center, double *Scatter)
{ /* multivariate Cauchy random generation to be called in R by 'rmCauchy' */
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  int info = 0, job = 1, n = *nobs, p = *nvar;

  GetRNGstate();
  chol_decomp(Scatter, p, p, job, &info);
  if (info)
    error("cholesky factorization in RNG_cauchy gave code %d", info);
  rmcauchy_std(y, n, p);
  mult_triangular_mat(1.0, Scatter, p, p, n, side, uplo, trans, diag, y, p);
  for (int i = 0; i < n; i++) {
    ax_plus_y(1., center, 1, y, 1, p);
    y += p;
  }
  PutRNGstate();
}

void
rmcauchy_std(double *y, int n, int p)
{ /* spherical Cauchy deviates */
  double tau, radial;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    tau = rgamma(.5, 2.);
    radial = R_pow(tau, -.5);
    scale_inplace(y, p, 1, radial);
    y += p;
  }
}

/* ========================================================================== *
 * multivariate contaminated normal random generation
 * ========================================================================== */

void
RNG_contaminated(double *y, int *nobs, int *nvar, double *center, double *Scatter, double *eps, double *vif)
{ /* multivariate contaminated normal random generation to be called in R 
   * by 'rmcnorm' */
  char *side = "L", *uplo = "U", *trans = "T", *diag = "N";
  int info = 0, job = 1, n = *nobs, p = *nvar;

  GetRNGstate();
  chol_decomp(Scatter, p, p, job, &info);
  if (info)
    error("cholesky factorization in RNG_contaminated gave code %d", info);
  rmcontaminated_std(y, *eps, *vif, n, p);
  mult_triangular_mat(1.0, Scatter, p, p, n, side, uplo, trans, diag, y, p);
  for (int i = 0; i < n; i++) {
    ax_plus_y(1., center, 1, y, 1, p);
    y += p;
  }
  PutRNGstate();
}

static void
rmcontaminated_std(double *y, double eps, double vif, int n, int p)
{ /* spherical contaminated normal deviates */
  double radial, unif;

  radial = 1. / vif;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++)
      y[j] = norm_rand();
    unif = unif_rand();
    if (unif < eps)
      scale_inplace(y, p, 1, radial);
    y += p;
  }
}
