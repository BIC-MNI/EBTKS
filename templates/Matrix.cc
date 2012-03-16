/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: Matrix.cc,v $
$Revision: 1.7 $
$Author: rotor $
$Date: 2010-05-18 23:01:21 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "Matrix.h"
#include <sys/stat.h>
#include "unistd.h"
#ifdef HAVE_MATLAB
#include "matlabSupport.h"
#endif
#include "miscTemplateFunc.h"

using namespace std;


#ifdef USE_COMPMAT
#ifdef __GNUC__
#pragma do_not_instantiate Mat<dcomplex>::array
#pragma do_not_instantiate Mat<dcomplex>::applyIndexFunction
#endif
#endif

#ifdef USE_FCOMPMAT
#ifdef __GNUC__
#pragma do_not_instantiate Mat<fcomplex>::array
#endif
#endif


/*********************************************************
template <class Type> const Mat<Type>::MATLAB = 0;
template <class Type> const Mat<Type>::RAW = 1;
template <class Type> const Mat<Type>::ASCII  = 2;
//template <class Type>  tells belong to template class 
*********************************************************/

#ifndef __GNUC__
template <class Type> unsigned Mat<Type>::_rangeErrorCount = 25;
template <class Type> Boolean  Mat<Type>::flushToDisk = FALSE;
#endif


#ifdef __GNUC__

/* For the GNU compiler, include the explicit specializations in this file,
 * so they can be seen during class instantiation.
 */
#include "MatrixSpec.cc"

#define _INSTANTIATE_MAT(Type) \
         template class Mat<Type>;                   \
         template Mat<Type>& pmultEquals(Mat<Type>&, const Mat<Type> &); \
         template Mat<Type>& pdivEquals(Mat<Type>&, const Mat<Type> &); \
         template Mat<Type> inv(const Mat<Type> &); \
         template Mat<Type> operator*<Type>(double, const Mat<Type> &); \
         template<> unsigned Mat<Type>::_rangeErrorCount = 25; \
         template<> Boolean  Mat<Type>::flushToDisk = FALSE;

_INSTANTIATE_MAT(int);
_INSTANTIATE_MAT(float);
_INSTANTIATE_MAT(double);
template class Zeros<double>;
template class Ones<double>;
template class Eye<double>;
template ostream & operator<<(ostream &, Mat<double> const &);
template SimpleArray<double> array(Mat<double> const &, double, double);
#ifdef USE_COMPMAT
_INSTANTIATE_MAT(dcomplex);
template Mat<dcomplex> fft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> fft(const Mat<dcomplex> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<double> &, unsigned, unsigned);
template Mat<dcomplex> ifft(const Mat<dcomplex> &, unsigned, unsigned);
#endif // USE_COMPMAT
#else
#ifdef USE_COMPMAT
template Mat<dcomplex>& pmultEquals(Mat<dcomplex>&, const Mat<dcomplex> &);
template Mat<dcomplex>& pdivEquals(Mat<dcomplex>&, const Mat<dcomplex> &);
#endif // USE_COMPMAT
template Mat<double>& pmultEquals(Mat<double>&, const Mat<double> &);
template Mat<double>& pdivEquals(Mat<double>&, const Mat<double> &);
#endif // __GNUC__ not defined
