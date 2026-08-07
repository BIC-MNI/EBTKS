/* Backend "cblas": the eight BLAS routines bundled LAPACK (clapack/dsytrf.c, dsytf2.c,
 * dlasyf.c, dsytrs.c) calls internally, forwarded to a named CBLAS instead of the bundled
 * clapack/*.c implementations. dsysv_ itself still comes from clapack/dsysv.c -- this
 * backend does not touch LAPACK, only the BLAS layer underneath it.
 *
 * The symbol names below are NOT dcopy_/dgemm_/etc. clapack/blaswrap.h -- included by every
 * bundled clapack/*.c file, unconditionally unless NO_BLAS_WRAP is defined, which nothing in
 * this tree does -- #defines dcopy_, dgemm_, dgemv_, dger_, dscal_, dswap_, dsyr_ and idamax_
 * to EBTKS_dcopy, EBTKS_dgemm, etc. before the bundled LAPACK/support files
 * (dsytrf.c/dsytf2.c/dlasyf.c/dsytrs.c) reference them, precisely so a real BLAS/LAPACK can
 * be linked into the same program without a multiple-definition clash on the plain names --
 * dsysv_ itself is deliberately left unrenamed, since that is the one symbol
 * legacy/N3/src/Splines/TBSpline.cc:648 calls from outside clapack/ entirely. This shim
 * therefore defines the EBTKS_-prefixed names directly; it does not include blaswrap.h.
 */

#include <assert.h>
#include <cblas.h>

typedef long int ebtks_integer;
typedef double   ebtks_doublereal;

static CBLAS_TRANSPOSE ebtks_cblas_trans(char c)
{
  switch(c)
    {
    case 'n': case 'N': return CblasNoTrans;
    case 't': case 'T': return CblasTrans;
    case 'c': case 'C': return CblasConjTrans;
    default:            assert(0 && "unrecognized transpose flag"); return CblasNoTrans;
    }
}

static CBLAS_UPLO ebtks_cblas_uplo(char c)
{
  switch(c)
    {
    case 'u': case 'U': return CblasUpper;
    case 'l': case 'L': return CblasLower;
    default:            assert(0 && "unrecognized uplo flag"); return CblasUpper;
    }
}

static int ebtks_narrow(ebtks_integer v)
{
  assert((ebtks_integer)(int) v == v && "value does not fit in int for this CBLAS call");
  return (int) v;
}

int EBTKS_dcopy(ebtks_integer *n, ebtks_doublereal *dx, ebtks_integer *incx,
                ebtks_doublereal *dy, ebtks_integer *incy)
{
  cblas_dcopy(ebtks_narrow(*n), dx, ebtks_narrow(*incx), dy, ebtks_narrow(*incy));
  return 0;
}

int EBTKS_dscal(ebtks_integer *n, ebtks_doublereal *da, ebtks_doublereal *dx,
                ebtks_integer *incx)
{
  cblas_dscal(ebtks_narrow(*n), *da, dx, ebtks_narrow(*incx));
  return 0;
}

int EBTKS_dswap(ebtks_integer *n, ebtks_doublereal *dx, ebtks_integer *incx,
                ebtks_doublereal *dy, ebtks_integer *incy)
{
  cblas_dswap(ebtks_narrow(*n), dx, ebtks_narrow(*incx), dy, ebtks_narrow(*incy));
  return 0;
}

/* clapack/idamax.c is 1-based (Fortran IDAMAX convention: 0 for n < 1, else the 1-based
 * position of the largest |element|); cblas_idamax returns a 0-based CBLAS_INDEX (size_t,
 * cblas.h). This is the single most likely silent defect in this backend if missed. */
ebtks_integer EBTKS_idamax(ebtks_integer *n, ebtks_doublereal *dx, ebtks_integer *incx)
{
  if(*n < 1 || *incx <= 0)
    return 0;
  return (ebtks_integer) cblas_idamax(ebtks_narrow(*n), dx, ebtks_narrow(*incx)) + 1;
}

int EBTKS_dgemv(char *trans, ebtks_integer *m, ebtks_integer *n, ebtks_doublereal *alpha,
                ebtks_doublereal *a, ebtks_integer *lda, ebtks_doublereal *x,
                ebtks_integer *incx, ebtks_doublereal *beta, ebtks_doublereal *y,
                ebtks_integer *incy)
{
  cblas_dgemv(CblasColMajor, ebtks_cblas_trans(*trans), ebtks_narrow(*m), ebtks_narrow(*n),
              *alpha, a, ebtks_narrow(*lda), x, ebtks_narrow(*incx), *beta, y,
              ebtks_narrow(*incy));
  return 0;
}

int EBTKS_dger(ebtks_integer *m, ebtks_integer *n, ebtks_doublereal *alpha,
               ebtks_doublereal *x, ebtks_integer *incx, ebtks_doublereal *y,
               ebtks_integer *incy, ebtks_doublereal *a, ebtks_integer *lda)
{
  cblas_dger(CblasColMajor, ebtks_narrow(*m), ebtks_narrow(*n), *alpha, x,
             ebtks_narrow(*incx), y, ebtks_narrow(*incy), a, ebtks_narrow(*lda));
  return 0;
}

int EBTKS_dsyr(char *uplo, ebtks_integer *n, ebtks_doublereal *alpha, ebtks_doublereal *x,
               ebtks_integer *incx, ebtks_doublereal *a, ebtks_integer *lda)
{
  cblas_dsyr(CblasColMajor, ebtks_cblas_uplo(*uplo), ebtks_narrow(*n), *alpha, x,
             ebtks_narrow(*incx), a, ebtks_narrow(*lda));
  return 0;
}

int EBTKS_dgemm(char *transa, char *transb, ebtks_integer *m, ebtks_integer *n,
                ebtks_integer *k, ebtks_doublereal *alpha, ebtks_doublereal *a,
                ebtks_integer *lda, ebtks_doublereal *b, ebtks_integer *ldb,
                ebtks_doublereal *beta, ebtks_doublereal *c__, ebtks_integer *ldc)
{
  cblas_dgemm(CblasColMajor, ebtks_cblas_trans(*transa), ebtks_cblas_trans(*transb),
              ebtks_narrow(*m), ebtks_narrow(*n), ebtks_narrow(*k), *alpha, a,
              ebtks_narrow(*lda), b, ebtks_narrow(*ldb), *beta, c__, ebtks_narrow(*ldc));
  return 0;
}
