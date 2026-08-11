/* Backend "lapacke": Cholesky solve over LAPACKE.
 *
 * Companion to shim/dsysv_lapacke.c.  TBSpline::fit assembles
 *
 *     A = lambda * nsamples * BendingEnergy + AtA
 *
 * a sum of two positive semidefinite matrices, hence positive semidefinite,
 * and positive definite whenever their null spaces meet only at zero.  dsysv
 * treats such a matrix as indefinite and pivots; the pivot sequence is chosen
 * by comparing floating-point magnitudes, so different LAPACK implementations
 * can choose differently and return solutions that differ by far more than
 * rounding once the matrix is ill-conditioned.  Cholesky does no pivoting, so
 * its factorization is determined by A alone.  It is also roughly half the
 * arithmetic and needs no workspace.
 *
 * Naming follows dsysv_lapacke.c and for the same reason: LAPACKE_dposv_work
 * calls dposv_ internally, so a dposv_ symbol defined here would be interposed
 * over the library's own and recurse.  Hence EBTKS_dposv.
 *
 * Integer width follows lapack_int, so ILP64 is handled here rather than at
 * the call site, which continues to use long int throughout.
 *
 * Called only from N3/src/Splines/TBSpline.cc, and only in the compilation
 * that defines N3_SPLINE_MODERN_SOLVE -- see N3/CMakeLists.txt.
 */

#include <lapacke.h>

typedef long int ebtks_integer;
typedef double   ebtks_doublereal;

/* Solve A x = b for symmetric positive definite A.  A is overwritten by its
 * Cholesky factor and b by the solution, matching dposv's contract.
 * info > 0 means the leading minor of that order is not positive definite,
 * i.e. A is not numerically SPD; the caller falls back to dsysv. */
int EBTKS_dposv(char *uplo, ebtks_integer *n, ebtks_integer *nrhs,
                ebtks_doublereal *a, ebtks_integer *lda,
                ebtks_doublereal *b, ebtks_integer *ldb, ebtks_integer *info)
{
  lapack_int info_ = LAPACKE_dposv_work(LAPACK_COL_MAJOR, *uplo,
                                        (lapack_int) *n, (lapack_int) *nrhs,
                                        a, (lapack_int) *lda,
                                        b, (lapack_int) *ldb);
  *info = (ebtks_integer) info_;
  return 0;
}
