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
$RCSfile: SimpleArray.cc,v $
$Revision: 1.9 $
$Author: bert $
$Date: 2004-12-08 17:05:18 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "EBTKS/SimpleArray.h"
#include "EBTKS/ValueMap.h"
#include <assert.h>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>

#ifdef HAVE_MATLAB
#include "matlabSupport.h"
#endif

#include "EBTKS/trivials.h"

namespace std //VF:fix for missing functions
{
  unsigned int abs(unsigned int v) { return v;}
  unsigned long abs(unsigned long v) { return v;}
  unsigned char abs(unsigned char v) { return v;}
};

using namespace std;

#ifdef USE_COMPMAT
#include "EBTKS/dcomplex.h"
#pragma do_not_instantiate SimpleArray<dcomplex>::log
#endif
#ifdef USE_FCOMPMAT
#include "EBTKS/fcomplex.h"
#pragma do_not_instantiate SimpleArray<fcomplex>::log
#endif

/*******************
 * SimpleArray class
 *******************/
// #ifndef __GNUC__
// template <class Type> unsigned SimpleArray<Type>::_rangeErrorCount = 25;
// #endif

// #ifdef __GNUC__
// #include "SimpleArraySpec.cc"
// 
// #define _INSTANTIATE_SIMPLEARRAY(Type)                    \
//          template class SimpleArray<Type>;                \
//          template class IndexStruct<Type>;                \
//          template ostream& operator << (ostream&, const SimpleArray<Type>&); \
//          template istream& operator >> (istream&, SimpleArray<Type>&); \
//          template SimpleArray<float> asFloatArray(const SimpleArray<Type>&); \
//          template SimpleArray<double> asDblArray(const SimpleArray<Type>&); \
//          template<> unsigned SimpleArray<Type>::_rangeErrorCount = 25;
// 
// #define _INSTANTIATE_SIMPLEARRAY_COMPLEX(Type)            \
//          template class SimpleArray<Type>;                \
//          template class IndexStruct<Type>;                \
//          template ostream& operator << (ostream&, const SimpleArray<Type>&); \
//          template istream& operator >> (istream&, SimpleArray<Type>&); \
//          template SimpleArray<double> asDblArray(const SimpleArray<Type>&); \
//          template<> unsigned SimpleArray<Type>::_rangeErrorCount = 25;
// 
// _INSTANTIATE_SIMPLEARRAY(char);
// _INSTANTIATE_SIMPLEARRAY(unsigned char);
// _INSTANTIATE_SIMPLEARRAY(short);
// _INSTANTIATE_SIMPLEARRAY(unsigned short);
// _INSTANTIATE_SIMPLEARRAY(int);
// _INSTANTIATE_SIMPLEARRAY(unsigned int);
// _INSTANTIATE_SIMPLEARRAY(float);
// _INSTANTIATE_SIMPLEARRAY(double);
// 
// 
// template SimpleArray<char> operator^(double, SimpleArray<char> const&);
// template SimpleArray<short> operator^(double, SimpleArray<short> const&);
// template SimpleArray<unsigned short> operator^(double, SimpleArray<unsigned short> const&);
// template SimpleArray<int> operator^(double, SimpleArray<int> const&);
// template SimpleArray<unsigned int> operator^(double, SimpleArray<unsigned int> const&);
// template SimpleArray<float> operator^(double, SimpleArray<float> const&);
// template SimpleArray<double> operator^(double, SimpleArray<double> const&);

// #ifdef USE_COMPMAT
// _INSTANTIATE_SIMPLEARRAY_COMPLEX(dcomplex);
// #endif // USE_COMPMAT
// 
// #ifdef USE_FCOMPMAT
// _INSTANTIATE_SIMPLEARRAY_COMPLEX(fcomplex);
// #endif // USE_FCOMPMAT

// #endif // __GNUC__

