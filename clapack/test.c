typedef long int _integer;
typedef float _real;
typedef double _doublereal;

extern int dsysv_(char *uplo, _integer *n, _integer *nrhs, _doublereal
        *a, _integer *lda, _integer *ipiv, _doublereal *b, _integer *ldb,
        _doublereal *work, _integer *lwork, _integer *info);

main()
{
  char uplo = 'u';
  _integer n = 10;
  _integer nrhs = 1;
  _integer lda = n;
  _integer ldb = n;
  _integer *ipiv = (_integer *) malloc(n * sizeof( _integer));
  _doublereal work;
  _integer lwork = 1;
  _doublereal *A;
  _doublereal *b;
  _integer info;

  dsysv_(&uplo, &n, &nrhs, (_doublereal *) A, &lda, ipiv, 
	 (_doublereal *) b, &ldb,
	 &work, &lwork, (_integer *)info);

  free(ipiv);

  return(b);
}

