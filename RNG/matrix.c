/* ID: matrix.c, last updated 2025-06-26, F.Osorio */

#include "RNG.h"

void
ax_plus_y(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y (AXPY operation) */

  /* quick return if possible */
  if (n <= 0 || incx <= 0 || incy <= 0)
    return;
  if (alpha == 0.0)
    return;

  if (incx == 1 && incy == 1) {
    /* code for increments equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      y[i] += alpha * x[i];

    for (int i = m; i + 7 < n; i += 8) {
      y[i] += alpha * x[i];
      y[i + 1] += alpha * x[i + 1];
      y[i + 2] += alpha * x[i + 2];
      y[i + 3] += alpha * x[i + 3];
      y[i + 4] += alpha * x[i + 4];
      y[i + 5] += alpha * x[i + 5];
      y[i + 6] += alpha * x[i + 6];
      y[i + 7] += alpha * x[i + 7];
    }
  } else {
    /* code for increments not equal to 1 */
    int ix = OFFSET(n, incx);
    int iy = OFFSET(n, incy);

    for (int i = 0; i < n; i++) {
      y[iy] += alpha * x[ix];
      ix += incx;
      iy += incy;
    }
  }
}

void
scale_inplace(double *x, int n, int inc, double alpha)
{ /* x <- alpha * x (x is overwritten) */

  /* quick return if possible */
  if (n <= 0 || inc <= 0)
    return;

  if (inc == 1) {
    /* code for increment equal to 1 */
    const int m = n % 8;

    for (int i = 0; i < m; i++)
      x[i] *= alpha;

    for (int i = m; i + 7 < n; i += 8) {
      x[i] *= alpha;
      x[i + 1] *= alpha;
      x[i + 2] *= alpha;
      x[i + 3] *= alpha;
      x[i + 4] *= alpha;
      x[i + 5] *= alpha;
      x[i + 6] *= alpha;
      x[i + 7] *= alpha;
    }
  } else {
    /* code for increment not equal to 1 */
    int ix = OFFSET(n, inc);

    for (int i = 0; i < n; i++) {
      x[ix] *= alpha;
      ix += inc;
    }
  }
}

void
mult_triangular_mat(double alpha, double *a, int lda, int nrow, int ncol, char *side, char *uplo, char *trans, char *diag, double *y, int ldy)
{ /* y <- alpha * op(a) %*% y, or y <- alpha * y %*% op(a),
   * with op(x) = x, or op(x) = t(x), and 'a' upper or lower triangular matrix */
  F77_CALL(dtrmm)(side, uplo, trans, diag, &nrow, &ncol, &alpha, a, &lda, y, &ldy FCONE FCONE FCONE FCONE);
}

void
chol_decomp(double *a, int lda, int p, int job, int *info)
{ /* cholesky factorization of a real symmetric positive definite matrix a.
   * the factorization has the form:
   * a <- l %*% t(l) (job = 0), or a <- t(u) %*% u (job = 1),
   * where u is an upper triangular matrix and l is lower triangular */
  char *uplo;

  uplo = (job) ? "U" : "L";
  F77_CALL(dpotrf)(uplo, &p, a, &lda, info FCONE);
}

/* DEBUG routine */

void
print_mat(double *x, int ldx, int nrow, int ncol, char *msg)
{ /* print matrix and message (used for printf debugging) */
  Rprintf( "%s\n", msg);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++)
      Rprintf( " %10.5g", x[i + j * ldx ]);
    Rprintf( "\n" );
  }
  Rprintf( "\n" );
}
