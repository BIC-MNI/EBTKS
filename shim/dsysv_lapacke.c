/* Backend "lapacke": dsysv_ implemented over LAPACKE, for a system LAPACK that ships the C
 * interface but whose Fortran-mangled symbols this shim prefers not to depend on directly.
 *
 * Signature matches legacy/N3/src/Splines/TBSpline.cc:70 exactly -- no header there changes.
 * That call site passes lwork = 1 with a one-element work array, not a workspace query
 * (:642-648), so LAPACKE_dsysv_work is used deliberately: LAPACKE_dsysv (no _work suffix)
 * ignores the caller's lwork and allocates its own optimal workspace internally, taking the
 * blocked path instead of the small unblocked one TBSpline.cc actually asks for.
 *
 * LAPACK_COL_MAJOR is mandatory: TBSpline.cc's A/B are laid out for the Fortran-ABI call
 * (column-major), and LAPACKE's row-major mode would transpose-copy them, changing the
 * arithmetic. Column-major is a pass-through.
 *
 * ipiv is allocated at lapack_int width (which follows the library this shim is built
 * against, so ILP64 is safe here even though the lapack backend has to declare it
 * unsupported) and copied back into the caller's long int[n]. TBSpline.cc never reads
 * ipiv after the call, but the copy costs nothing and removes the dependency on that
 * observation. info is likewise narrowed/widened through the long int* the caller passed.
 */

#include <lapacke.h>
#include <stdlib.h>

typedef long int ebtks_integer;
typedef double   ebtks_doublereal;

int dsysv_(char *uplo, ebtks_integer *n, ebtks_integer *nrhs, ebtks_doublereal *a,
           ebtks_integer *lda, ebtks_integer *ipiv, ebtks_doublereal *b,
           ebtks_integer *ldb, ebtks_doublereal *work, ebtks_integer *lwork,
           ebtks_integer *info)
{
  lapack_int n_ = (lapack_int) *n;
  lapack_int nrhs_ = (lapack_int) *nrhs;
  lapack_int lda_ = (lapack_int) *lda;
  lapack_int ldb_ = (lapack_int) *ldb;
  lapack_int lwork_ = (lapack_int) *lwork;
  lapack_int *ipiv_ = NULL;
  lapack_int info_;
  ebtks_integer i;

  if(n_ > 0)
    ipiv_ = (lapack_int *) malloc(sizeof(lapack_int) * (size_t) n_);

  info_ = LAPACKE_dsysv_work(LAPACK_COL_MAJOR, *uplo, n_, nrhs_, a, lda_, ipiv_,
                              b, ldb_, work, lwork_);

  for(i = 0; i < *n; i++)
    ipiv[i] = (ebtks_integer) ipiv_[i];

  free(ipiv_);
  *info = (ebtks_integer) info_;

  /* f2c-generated dsysv_ always returns 0; TBSpline.cc consults *info, not this value. */
  return 0;
}
