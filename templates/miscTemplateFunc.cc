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
$RCSfile: miscTemplateFunc.cc,v $
$Revision: 1.2 $
$Author: jason $
$Date: 2002-03-20 21:42:47 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "miscTemplateFunc.h"
#include "Histogram.h"
#include "dcomplex.h"
#include "fcomplex.h"

double
asDouble(dcomplex value) 
{
  return sqrt(norm(value));
}

double
asDouble(fcomplex value) 
{
  return sqrt(norm(value));
}

template <class Type>
double
asDouble(Type value) 
{
  return double(value);
}

template <class Type>
Histogram&
add(Histogram& hist, const SimpleArray<Type>& array)
{
  if (array.size()) {
    array.resetIterator();
    for (unsigned i = array.size(); i; i--)
      hist.add(array++);
  }

  return hist;
}

template <class Type>
Type nextPowerOf2(Type value)
{
  double x;

  if (value >= 0) {
    x = 1;
    while (x < value)
      x *= 2;
  }
  else {
    x = -1;
    while (x > value)
      x *= 2;
  }

  return Type(x);
}

template <class Type>
int isPowerOf2(Type value)
{
  if (!value)
    return 0;

  if (value < 0)
    value = -value;

  Type x = 1;
  while (x < value)
    x *= 2;

  return x == value;
}

#ifdef __GNUC__
template Histogram& add(Histogram& hist, const SimpleArray<float>& array);

#define _INSTANTIATE_MISCFUNC(Type) \
   template double asDouble(Type);        \
   template void swap(Type &, Type &);  \
   template Type min(const Type &, const Type &); \
   template Type max(const Type &, const Type &); \
   template Type min(const Type &, const Type &, const Type &); \
   template Type max(const Type &, const Type &, const Type &); \
   template Type clamp(const Type &, const Type &, const Type &); \
   template Type nextPowerOf2(Type); \
   template int isPowerOf2(Type);

_INSTANTIATE_MISCFUNC(char);
_INSTANTIATE_MISCFUNC(unsigned char);
_INSTANTIATE_MISCFUNC(short);
_INSTANTIATE_MISCFUNC(int);
_INSTANTIATE_MISCFUNC(unsigned int);
_INSTANTIATE_MISCFUNC(float);
_INSTANTIATE_MISCFUNC(double);
template void swap(dcomplex &, dcomplex &);  
#endif
