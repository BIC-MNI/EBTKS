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
$RCSfile: ValueMap.h,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:26 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _VALUE_MAP_H
#define _VALUE_MAP_H

#include <iostream.h>
#include "SimpleArray.h"

/*************************
 * Abstract ValueMap class
 *************************/

class ValueMap {
public:
  // Invert the value map
  virtual ValueMap& inv() = 0;

  // Map concatenation
  virtual ValueMap& concat(const ValueMap& map) = 0;
  virtual ValueMap& operator () (const ValueMap& map);

  // Evaluate map
  virtual double operator () (double sourceValue) const = 0;
  // Reverse evaluate map
  virtual double reverse(double destValue) const = 0;

  // I/O
  virtual ostream& print(ostream&) const = 0;
};

ostream& operator << (ostream& os, const ValueMap& map);

/******************
 * Linear map class
 ******************/

class LinearMap : public ValueMap {
  double _factor;
  double _offset;

public:
  LinearMap(double factor = 1.0, double offset = 0.0) { 
    _factor = factor; _offset = offset; }
  LinearMap(const LinearMap& map) { 
    _factor = map._factor; _offset = map._offset; }
  LinearMap(double sourceMin, double sourceMax, double destMin, double destMax) {
    _factor = (destMax - destMin)/(sourceMax - sourceMin);
    _offset = destMin - _factor*sourceMin; }
  virtual ~LinearMap() { _factor = 1.0; _offset = 0.0; }

  double factor() const { return _factor; }
  double offset() const { return _offset; }

  double& factor() { return _factor; }
  double& offset() { return _offset; }

  double factor(double f) { _factor = f; return _factor; }
  double offset(double o) { _offset = o; return _offset; }

  LinearMap& operator = (const LinearMap& map) { 
    _factor = map._factor; _offset = map._offset; return *this; }

  ValueMap& inv();
  ValueMap& concat(const ValueMap& map);
  ValueMap& concat(const LinearMap& map);

  LinearMap& operator () (double factor, double offset);
  LinearMap& operator () (double sourceMin, double sourceMax, double destMin,
			  double destMax);
  
  double operator () (double sourceValue) const { return _offset+_factor*sourceValue;}
  double reverse(double destValue) const        { return (destValue-_offset)/_factor;}
  
  ostream& print(ostream& os) const;
};

/********************
 * Lookup table class
 ********************/

template <class Type>
class LUT : public ValueMap {
  SimpleArray<Type> _source;
  SimpleArray<Type> _dest;

public:
  LUT(unsigned length = 0); // Allocated size; actual size is zero
  LUT(const SimpleArray<Type>& source, const SimpleArray<Type>& dest);
  LUT(const LUT& map);
  virtual ~LUT();

  LUT& add(Type source, Type dest); // Add one map entry

  LUT& operator = (const LUT& map);

  ValueMap& inv();
  ValueMap& concat(const ValueMap& map);
  ValueMap& concat(const LUT& map);

  double operator () (double sourceValue) const;
  double reverse(double destValue) const;
  
  ostream& print(ostream&) const;

private:
  void _sort();
};

/**************************
 * Mapping of other objects
 **************************/

// Map all elements through a single map
template <class Type>
SimpleArray<Type>& map(SimpleArray<Type>& array, const ValueMap&);
// Removed because of DCC resolving problems
//template <class Type>
//SimpleArray<Type>  map(const SimpleArray<Type>& array, const ValueMap& valueMap) {
//  return mapConst(array, valueMap); }
template <class Type>
SimpleArray<Type>  mapConst(const SimpleArray<Type>& array, const ValueMap& valueMap);

// Map each element through a linear map
template <class Type>
SimpleArray<Type>& map(SimpleArray<Type>& array, const Array<LinearMap>&);
// Removed because of DCC resolving problems
//template <class Type>
//SimpleArray<Type>  map(const SimpleArray<Type>& array, const Array<LinearMap>& maps) {
//  return mapConst(array, maps); }
template <class Type>
SimpleArray<Type>  mapConst(const SimpleArray<Type>& array,const Array<LinearMap>& maps);

#endif

