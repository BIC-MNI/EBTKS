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
$RCSfile: miscTemplateFunc.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:26 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MISC_TEMPLATE_FUNC_H
#define _MISC_TEMPLATE_FUNC_H

#include <stdlib.h>

// Forward declarations
class Histogram;
template <class Type> class SimpleArray;

// Double cast which works for complex data
template <class Type>
double asDouble(Type value);

// Value swapping
template <class Type>
inline void swap(Type& x, Type& y) { Type tmp(x); x = y; y = tmp; }

// Power functions
inline double intPower(double x, int y) {
  if (!y) return(1.0);
  if (!x) return(0.0);
  if (x == 1.0) return(1.0);

  register double result = x;
  for (register unsigned i = abs(y) - 1; i; i--) 
    result *= x;

  return (y < 0) ? 1.0/result : (double) result;
}

template <class Type>
Type nextPowerOf2(Type value);
template <class Type>
int isPowerOf2(Type value);

// Two argument min/max
template <class Type>
inline Type min(const Type& x, const Type& y) { return (x < y) ? x : y; }
template <class Type>
inline Type max(const Type& x, const Type& y) { return (x > y) ? x : y; }

// Three argument min/max
template <class Type>
inline Type min(const Type& x, const Type& y, const Type& z) { 
  Type value = (x < y) ? x : y; return (value < z) ? value : z; }
template <class Type>
inline Type max(const Type& x, const Type& y, const Type& z) { 
  Type value = (x > y) ? x : y; return (value > z) ? value : z; }

// Clamping
template <class Type>
inline Type clamp(const Type& value, const Type& a, const Type& b) {
  if (b < a) {
    if (value > a)
      return a;
    if (value < b)
      return b;
  }
  else {
    if (value < a)
      return a;
    if (value > b)
      return b;
  }

  return value;
}

template <class Type>
Histogram& add(Histogram& hist, const SimpleArray<Type>& array);

#endif
