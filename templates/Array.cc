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
$RCSfile: Array.cc,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2004-12-08 16:43:44 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "EBTKS/Array.h"
#include "EBTKS/dcomplex.h"
#include <assert.h>
#include <string.h>
#include <iostream>		// (bert) changed from iostream.h
#include <math.h>
//#include <string>		// (bert) changed from string.h
using namespace std;		// (bert) added

/******************
 * Array base class
 ******************/

// #ifdef __GNUC__
// #define _INSTANTIATE_ARRAY(Type)                       \
//          template class Array<Type>;                   \
//          template<> unsigned int Array<Type>::_arrayCtr = 0;          \
//          template<> Boolean  Array<Type>::_debug = FALSE;         \
//          template<> unsigned int Array<Type>::_rangeErrorCount = 25;  \
//          template<> unsigned size(const Array<Type> &);
// 
// _INSTANTIATE_ARRAY(char);
// _INSTANTIATE_ARRAY(unsigned char);
// _INSTANTIATE_ARRAY(short);
// _INSTANTIATE_ARRAY(int);
// _INSTANTIATE_ARRAY(unsigned int);
// _INSTANTIATE_ARRAY(unsigned short);
// _INSTANTIATE_ARRAY(float);
// _INSTANTIATE_ARRAY(double);
// _INSTANTIATE_ARRAY(dcomplex);
// _INSTANTIATE_ARRAY(fcomplex);
// #include "EBTKS/Path.h"
// _INSTANTIATE_ARRAY(Path);
// #include "EBTKS/ValueMap.h"
// _INSTANTIATE_ARRAY(LinearMap);
// #include "EBTKS/SimpleArray.h"
// _INSTANTIATE_ARRAY(SimpleArray<char> );
// _INSTANTIATE_ARRAY(SimpleArray<unsigned char> );
// _INSTANTIATE_ARRAY(SimpleArray<short> );
// _INSTANTIATE_ARRAY(SimpleArray<int> );
// _INSTANTIATE_ARRAY(SimpleArray<unsigned int> );
// _INSTANTIATE_ARRAY(SimpleArray<unsigned short> );
// _INSTANTIATE_ARRAY(SimpleArray<float> );
// _INSTANTIATE_ARRAY(SimpleArray<double> );
// _INSTANTIATE_ARRAY(Array< SimpleArray<double> >);
// #endif


// #ifndef __GNUC__
// template <class Type> unsigned Array<Type>::_arrayCtr = 0;
// template <class Type> Boolean  Array<Type>::_debug = FALSE;
// template <class Type> unsigned Array<Type>::_rangeErrorCount = 25;
// #endif


