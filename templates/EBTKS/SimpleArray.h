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
$RCSfile: SimpleArray.h,v $
$Revision: 1.5 $
$Author: bert $
$Date: 2004-12-08 17:05:18 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef SIMPLE_ARRAY_H
#define SIMPLE_ARRAY_H

#include <math.h>

#include "EBTKS/trivials.h"
#include "EBTKS/Array.h"

/********************************************************************
 * SimpleArray class
 *
 * This class is derived from Array. It extends the functionality of
 * Array with various mathematical operations; as a result, these
 * operations must be defined for the objects stored in the array.
 * SimpleArray is therefore most suitable for simple types.
 ********************************************************************/

class ValueMap;
template <class Type> class SimpleArray;

typedef SimpleArray<char>          BoolArray;
typedef SimpleArray<int>           IntArray;
typedef SimpleArray<unsigned>      UnsignedArray;
typedef SimpleArray<double>        DblArray;

template <class Type>
struct IndexStruct {
  Type     element;
  unsigned index;

  static int compareAscending(const void *struct1, const void *struct2);
  static int compareDescending(const void *struct1, const void *struct2);
};

template <class Type>
class SimpleArray : public Array<Type> {
public:
  static int     compareAscending(const void *, const void *);
  static int     compareDescending(const void *, const void *);
//   static unsigned     _rangeErrorCount;

// Manager functions
  SimpleArray(unsigned sz = DEFAULT_SIZE) : Array<Type>(sz) {}
  SimpleArray(Type value, unsigned sz) : Array<Type>(value, sz) {}
  SimpleArray(const Type *init, unsigned nElements) : Array<Type>(init, nElements) {}
  SimpleArray(const Array<Type>& array) : Array<Type>(array) {}
  SimpleArray(const SimpleArray<Type>& array) : Array<Type>(array) {}
  SimpleArray(Type minVal, Type step, Type maxVal);
  virtual ~SimpleArray() {}

  // Binary I/O functions
  std::ostream&         write(std::ostream& os) const      { saveBinary(os); return os; }
  virtual std::ostream& saveBinary(std::ostream&, unsigned n = 0, unsigned start = 0) const;
  std::istream&         read(std::istream& is)             { loadBinary(is); return is; }
  virtual std::istream& loadBinary(std::istream&, unsigned n = 0, unsigned start = 0);

  // ASCII I/O functions
  std::ostream&         print(std::ostream& os) const     { saveAscii(os); return os; }
  virtual std::ostream& saveAscii(std::ostream&, unsigned n = 0, unsigned start = 0) const;
  std::istream&         scan(std::istream& is)            { loadAscii(is); return is; }
  virtual std::istream& loadAscii(std::istream&, unsigned n = 0, unsigned start = 0);

  // Matlab I/O functions
#ifdef HAVE_MATLAB
  virtual Boolean saveMatlab(const char *fileName, const char *varName = "A",
			     const char *option = "u") const;
#endif

// Get functions
  SimpleArray operator () (unsigned nElements) const;   // Sub-array (first nElements)
  SimpleArray operator () (unsigned start, unsigned end) const; // Sub-array
  SimpleArray operator () (const BoolArray& boolArray) const;
  // Return an array containing only the elements for which boolArray is TRUE
  SimpleArray operator () (const UnsignedArray& indices) const;
  // Return an array containing only the elements with indices listed in <indices>

  Boolean  contains(Type value) const;
  Boolean  contains(Type value, unsigned start, unsigned end) const;
  Boolean  containsOnly(Type value) const;
  Boolean  containsOnly(Type value, unsigned start, unsigned end) const;
  // Total number of occurrences of object
  unsigned occurrencesOf(Type value) const {
    return occurrencesOf(value, 0, this->_size - 1); }
  // Total number of occurrences of object present in a range
  unsigned occurrencesOf(Type, unsigned start, unsigned end) const;
  // Index of first (dir > 0) or last (dir < 0) occurrence of value
  int    indexOf(Type, int dir, unsigned start) const;
  int    indexOf(Type value, int dir = 1) const {
    return indexOf(value, dir, (dir > 0) ? 0 : this->_size - 1); }
  int    indexOf(Type value, unsigned start) const {
    return indexOf(value, 1, start); }

  void  removeAll(Type value);
  // Remove all elements equal to <value>
  void  removeAllIn(Type floor, Type ceil, unsigned *N = 0);
  // Remove all elements in [floor, ceil]
  void  removeAllNotIn(Type floor, Type ceil, unsigned *nBelow=0, unsigned *nAbove=0);
  // Remove all elements not in [floor, ceil]
  UnsignedArray indicesOf(Type value) const;
  // Return an array containing the (unique) elements common with <array>
  SimpleArray common(const SimpleArray<Type>& array) const;

  // Removes all infinite and NaN values
  SimpleArray& prune();
  // Fills array with uniformly distributed numbers
  SimpleArray& randuniform(double min = 0, double max = 1);
  // Fills array with normally distributed numbers
  SimpleArray& randnormal(double mean = 0, double std = 1);

  // Quicksort all elements
  virtual void qsort() { qsortAscending(); }
  virtual void qsort(int (*compare) (const void *, const void *)) {
    ::qsort(this->_contents, this->_size, sizeof(Type), compare); }
  virtual void qsortAscending() {
    ::qsort(this->_contents, this->_size, sizeof(Type), compareAscending); }
  virtual void qsortDescending() {
    ::qsort(this->_contents, this->_size, sizeof(Type), compareDescending); }
  virtual SimpleArray<unsigned> qsortIndexAscending() const;
  virtual SimpleArray<unsigned> qsortIndexDescending() const;

  Type   min(unsigned *index = 0) const;
  Type   max(unsigned *index = 0) const;
  void   extrema(Type *min, Type *max) const;
  Type   range(unsigned *minIndex = 0, unsigned *maxIndex = 0) const;
  double sum() const;
  double sum2() const;
  double mean() const { return sum()/double(this->_size); }
  double prod() const;
  double prod2() const;
  double var() const;
  double std() const { return ::sqrt(var()); }

  // Median functions. median() and medianVolatile() use the algorithm described
  // in Sections 8.1, 8.3, and 10.2 in Cormen, Leierson, and Rivest,
  // "Introduction to Algorithms" (aka the Big White Book).
  // The low median is returned in case the array has an even # elements.
  virtual Type median() const;
  virtual Type medianVolatile(); // Does not create a copy of the array, but changes it

  virtual Type mode(const Type binWidth) const;

  // Cumulative sum/prod2uct
  DblArray cumSum() const;
  DblArray cumProd() const;

  // Impose ceiling/floor to values
  void ceil(Type ceil);
  void floor(Type floor);

  // Boolean operations
  // Elementwise inequality/equality
  Boolean operator != (const SimpleArray&) const;
  Boolean operator == (const SimpleArray& array) const {
    return !((*this) != array); }

  BoolArray operator && (const SimpleArray<Type>& array) const;
  BoolArray operator || (const SimpleArray<Type>& array) const;

  BoolArray operator == (double value) const;
  BoolArray operator != (double value) const;
  BoolArray operator >= (double value) const;
  BoolArray operator >  (double value) const;
  BoolArray operator <= (double value) const;
  BoolArray operator <  (double value) const;

  BoolArray operator >= (const SimpleArray<Type>&) const;
  BoolArray operator >  (const SimpleArray<Type>&) const;
  BoolArray operator <= (const SimpleArray<Type>&) const;
  BoolArray operator <  (const SimpleArray<Type>&) const;

// Arithmetic operations
  SimpleArray  operator - () const {
    SimpleArray<Type> result(Type(0), this->_size); return result -= *this; }

  SimpleArray& operator += (Type);
  SimpleArray  operator +  (Type value) const {
    SimpleArray<Type> result(*this); return result += value; }
  SimpleArray& operator += (const SimpleArray&);
  SimpleArray  operator +  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result += array; }

  SimpleArray& operator -= (Type);
  SimpleArray  operator -  (Type value) const {
    SimpleArray<Type> result(*this); return result -= value; }
  SimpleArray& operator -= (const SimpleArray&);
  SimpleArray  operator -  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result -= array; }

  SimpleArray& operator *= (double);
  SimpleArray  operator *  (double value) const {
    SimpleArray<Type> result(*this); return result *= value; }
  SimpleArray& operator *= (const SimpleArray&);
  SimpleArray  operator *  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result *= array; }

  SimpleArray& operator /= (double value) {
    return (*this) *= 1.0/value; }
  SimpleArray  operator /  (double value) const {
    SimpleArray<Type> result(*this); return result /= value; }
  SimpleArray& operator /= (const SimpleArray&);
  SimpleArray  operator /  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result /= array; }

  SimpleArray abs() const;                // Absolute value
  SimpleArray round(unsigned n = 0) const;// Round to n decimal places
  SimpleArray sqr() const;                   // Square of all elements
  SimpleArray sqrt() const;                  // Square root of all elements
  SimpleArray operator ^ (int power) const;  // Raise all elements to <power>
    // Note: all elements are cast to int first. Should be generalized!!
  SimpleArray ln() const;                    // Natural logarithm of all elements
  SimpleArray log() const;                   // Base-10 logarithm of all elements
  SimpleArray exp() const;                   // e^(all elements)
  SimpleArray exp10() const;                 // 10^(all elements)

  SimpleArray sample(unsigned maxN) const;
  SimpleArray applyElementWise(Type (*function) (Type)) const;
  SimpleArray map(const ValueMap& map) const;

protected:
  // Median support functions
  Type _randomizedSelect(int p, int r, int i);
  int  _randomizedPartition(int p, int r);
  int  _partition(int p, int r);
};

// Various non-member template functions
//
template <class Type>
unsigned size(const SimpleArray<Type>& array) { return array.size(); }

template <class Type>
Type min(const SimpleArray<Type>& array) { return array.min(); }

template <class Type>
Type max(const SimpleArray<Type>& array) { return array.max(); }

template <class Type>
Type range(const SimpleArray<Type>& array, unsigned *minIndex = 0,
	   unsigned *maxIndex = 0)
{
  return array.range(minIndex, maxIndex);
}

template <class Type>
double sum(const SimpleArray<Type>& array) { return array.sum(); }

template <class Type>
double sum2(const SimpleArray<Type>& array) { return array.sum2(); }

template <class Type>
double mean(const SimpleArray<Type>& array) { return array.mean(); }

template <class Type>
double prod(const SimpleArray<Type>& array) { return array.prod(); }

template <class Type>
double prod2(const SimpleArray<Type>& array) { return array.prod2(); }

template <class Type>
double var(const SimpleArray<Type>& array) { return array.var(); }

template <class Type>
double stdev(const SimpleArray<Type>& array) { return array.std(); }

template <class Type>
Type median(const SimpleArray<Type>& array) { return array.median(); }

template <class Type>
Type medianVolatile(SimpleArray<Type>& array) { return array.medianVolatile(); }

template <class Type>
Type mode(const SimpleArray<Type>& array, Type binWidth) { return array.mode(binWidth); }

template <class Type>
SimpleArray<Type> abs(const SimpleArray<Type>& array) { return array.abs(); }

template <class Type>
SimpleArray<Type> sqr(const SimpleArray<Type>& array) { return array.sqr(); }

template <class Type>
SimpleArray<Type> sqrt(const SimpleArray<Type>& array) { return array.sqrt(); }

template <class Type>
SimpleArray<Type> ln(const SimpleArray<Type>& array) { return array.ln(); }

template <class Type>
SimpleArray<Type> log(const SimpleArray<Type>& array) { return array.log(); }

template <class Type>
SimpleArray<Type> max(const SimpleArray<Type>& array1,
		      const SimpleArray<Type>& array2);

template <class Type>
SimpleArray<Type> min(const SimpleArray<Type>& array1,
		      const SimpleArray<Type>& array2);

template <class Type>
SimpleArray<double> cumSum(const SimpleArray<Type>& array)
{
  return array.cumSum();
}

template <class Type>
SimpleArray<double> cumProd(const SimpleArray<Type>& array)
{
  return array.cumProd();
}

template <class Type>
SimpleArray<Type> exp(const SimpleArray<Type>& array)
{
  return ::exp(1.0) ^ array;
}

template <class Type>
SimpleArray<Type> exp10(const SimpleArray<Type>& array) { return 10^array; }

template <class Type>
SimpleArray<Type> operator ^ (double base, const SimpleArray<Type>& array);

// I/O
template <class T> std::ostream& operator << (std::ostream& os, const SimpleArray<T>& A)
{
  return A.saveAscii(os);
}

template <class T> std::istream& operator >> (std::istream& is, SimpleArray<T>& A)
{
  return A.loadAscii(is);
}

// Type conversions
template <class Type>
SimpleArray<int> asIntArray(const SimpleArray<Type>&);

template <class Type>
SimpleArray<float> asFloatArray(const SimpleArray<Type>&);

template <class Type>
SimpleArray<double> asDblArray(const SimpleArray<Type>&);

#ifdef USE_COMPMAT
template <>
SimpleArray<dcomplex> SimpleArray<dcomplex>::round(unsigned n) const;

template <>
SimpleArray<dcomplex>& SimpleArray<dcomplex>::prune();

template <>
SimpleArray<dcomplex> operator ^ (double, const SimpleArray<dcomplex>&);
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template <>
SimpleArray<fcomplex> SimpleArray<fcomplex>::round(unsigned n) const;

template <>
SimpleArray<fcomplex>& SimpleArray<fcomplex>::prune();

template <>
SimpleArray<fcomplex> operator ^ (double, const SimpleArray<fcomplex>&);
#endif /* USE_FCOMPMAT */

template <class Type>
SimpleArray<Type> round(const SimpleArray<Type>& array, unsigned n = 0)
{
  return array.round(n);
}

template <class Type>
void prune(SimpleArray<Type>& array)
{
  array.prune();
}


template <class Type>
int
IndexStruct<Type>::compareAscending(const void *struct1, const void *struct2) {
  const Type& element1 = ((IndexStruct<Type> *) struct1)->element;
  const Type& element2 = ((IndexStruct<Type> *) struct2)->element;
  if (element1 > element2)
    return 1;
  if (element1 < element2)
    return -1;
  return 0;
}

template <class Type>
int
IndexStruct<Type>::compareDescending(const void *struct1, const void *struct2) {
  const Type& element1 = ((IndexStruct<Type> *) struct1)->element;
  const Type& element2 = ((IndexStruct<Type> *) struct2)->element;
  if (element1 < element2)
    return 1;
  if (element1 > element2)
    return -1;
  return 0;
}

template <class Type>
int
SimpleArray<Type>::compareAscending(const void *ptr1, const void *ptr2)
{
  const Type& element1 = *(Type *) ptr1;
  const Type& element2 = *(Type *) ptr2;
  if (element1 > element2)
    return 1;
  if (element1 < element2)
    return -1;
  return 0;
}

template <class Type>
int
SimpleArray<Type>::compareDescending(const void *ptr1, const void *ptr2)
{
  const Type& element1 = *(Type *) ptr1;
  const Type& element2 = *(Type *) ptr2;
  if (element1 < element2)
    return 1;
  if (element1 > element2)
    return -1;
  return 0;
}

template <class Type>
SimpleArray<Type>::SimpleArray(Type minVal, Type step, Type maxVal)
: Array<Type>(unsigned(std::abs((double)(maxVal - minVal)/step)) + 1)
{
  register Type *elementPtr = this->_contents;
  Type value = minVal;
  for (register unsigned i = this->_size; i; i--) {
    *elementPtr++ = value;
    value = Type(step + value);
  }
}

//
// Binary I/O functions
//

template <class Type>
std::ostream&
SimpleArray<Type>::saveBinary(std::ostream& os, unsigned n, unsigned start) const
{
  if (start >= this->_size) {
    if (this->_size)
//       if (_rangeErrorCount) {
//         _rangeErrorCount--;
        std::cerr << "SimpleArray::saveBinary: start out of range" << std::endl;
//       }

    return os;
  }

  if (!n)
    n = this->_size - start;
  else if (start + n > this->_size) {
    n = this->_size - start;
//     if (_rangeErrorCount) {
//       _rangeErrorCount--;
      std::cerr << "SimpleArray::saveBinary: n too large; truncated" << std::endl;
//     }
  }

  os.write((char *) this->_contents + start, n*sizeof(Type));

  return os;
}

template <class Type>
std::istream&
SimpleArray<Type>::loadBinary(std::istream& is, unsigned n, unsigned start)
{
  if (!n)
    n = this->_size;

  this->newSize(start + n);

  if (this->_size)
    is.read((char *) this->_contents + start, n*sizeof(Type));

  return is;
}

//
// Ascii I/O functions
//

template <class Type>
std::ostream&
SimpleArray<Type>::saveAscii(std::ostream& os, unsigned n, unsigned start) const
{
  if (start >= this->_size) {
    if (this->_size)
//       if (_rangeErrorCount) {
//         _rangeErrorCount--;
        std::cerr << "SimpleArray::saveAscii: start out of range" << std::endl;
//       }
    return os;
  }

  if (!n)
    n = this->_size - start;
  else if (start + n > this->_size) {
    n = this->_size - start;
//     if (_rangeErrorCount) {
//       _rangeErrorCount--;
      std::cerr << "SimpleArray::saveAscii: n too large; truncated" << std::endl;
//     }
  }

  this->resetIterator(start);
  for (unsigned i = n; i && os; i--) {
    os << (*this)++;
    if (i > 1)
      os << " ";
  }

  return os;
}

template <class Type>
std::istream&
SimpleArray<Type>::loadAscii(std::istream& is, unsigned n, unsigned start)
{
  if (!n)
    n = this->_size;

  this->newSize(start + n);

  this->resetIterator(start);
  for (unsigned i = n; i && is; i--)
    is >> (*this)++;

  return is;
}

#ifdef HAVE_MATLAB
Boolean
SimpleArray<double>::saveMatlab(const char *fileName, const char *varName,
        const char *option) const
{
  return ::saveMatlab(fileName, varName, option, _size, 1, this->_contents);
}

#ifdef USE_COMPMAT
Boolean
SimpleArray<dcomplex>::saveMatlab(const char *fileName, const char *varName,
         const char *option) const
{
  double *real = 0;
  double *imag = 0;

  if (_size) {
    real = new double[_size];
    if (!real) {
      std::cerr << "Couldn't allocate temporary real array for MATLAB conversion" << std::endl;
      return FALSE;
    }

    imag = new double[_size];
    if (!imag) {
      std::cerr << "Couldn't allocate temporary imag array for MATLAB conversion" << std::endl;
      return FALSE;
    }

    double        *realPtr = real;
    double        *imagPtr = imag;
    const dcomplex *contentsPtr = this->_contents;

    for (unsigned i = _size; i; i--) {
      *realPtr++ = ::real(*contentsPtr);
      *imagPtr++ = ::imag(*contentsPtr);
      contentsPtr++;
    }
  }

  Boolean status = ::saveMatlab(fileName, varName, option, _size, 1, real, imag);

  if (real)
    delete [] real;
  if (imag)
    delete [] imag;

  return status;
}
#endif

template <class Type>
Boolean
SimpleArray<Type>::saveMatlab(const char *fileName, const char *varName,
            const char *option) const
{
  double *temp = 0;
  if (_size) {
    temp = new double[_size];
    if (!temp) {
      std::cerr << "Couldn't allocate temporary double array for MATLAB conversion" << std::endl;
      return FALSE;
    }

    double     *tempPtr     = temp;
    const Type *contentsPtr = this->_contents;
    for (unsigned i = _size; i; i--)
      *tempPtr++ = asDouble(*contentsPtr++);
  }

  Boolean status = ::saveMatlab(fileName, varName, option, _size, 1, temp);

  if (temp)
    delete [] temp;

  return status;
}
#endif

//
// Get functions
//

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::operator () (unsigned nElements) const
{
  if (nElements > this->_size) {
    std::cerr << "Warning! Array::operator(" << nElements
   << ") called with on array of size " << this->_size << ". Value truncated!" << std::endl;
    nElements = this->_size;
  }

  SimpleArray<Type> subArray(nElements);

  Type *sourcePtr = this->_contents;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = nElements; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::operator () (unsigned start, unsigned end) const
{
  unsigned n = end - start + 1;
  if (start + n > this->_size) {
    std::cerr << "Warning! Array::operator(" << start << ", " << end
   << ") called with on array of size " << this->_size << ". Truncated!" << std::endl;
    n = this->_size - start;
  }

  SimpleArray<Type> subArray(n);

  Type *sourcePtr = this->_contents + start;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = n; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::operator () (const BoolArray& boolArray) const
{
  unsigned i;
  unsigned size = MIN(this->_size, boolArray.size());
  unsigned newSize = 0;
  const Boolean *boolPtr = boolArray.contents();
  for (i = size; i; i--)
    if (*boolPtr++)
      newSize++;

  SimpleArray<Type> array(newSize);
  boolPtr = boolArray.contents();
  Type    *sourcePtr = this->_contents;
  Type    *destPtr   = array._contents;
  for (i = size; i; i--) {
    if (*boolPtr++)
      *destPtr++ = *sourcePtr;
    sourcePtr++;
  }

  return array;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::operator () (const UnsignedArray& indices) const
{
  unsigned trueN = 0;
  unsigned N     = indices.size();
  SimpleArray<Type> array(N);
  Type     *destPtr = array.contents();
  const unsigned *index   = indices.contents();
  for (unsigned i = N; i; i--, index++) {
    if (*index >= this->_size) {
//       if (_rangeErrorCount) {
//   _rangeErrorCount--;
        std::cerr << "Warning! SimpleArray::operator(): index " << *index
       << "out of range!" << std::endl;
//       }
    }
    else {
      *destPtr++ = this->_contents[*index];
      trueN++;
    }
  }

  array.newSize(trueN);

  return array;
}

template <class Type>
Boolean
SimpleArray<Type>::contains(Type value) const
{
  register Type *contentsPtr = this->_contents;
  for (register unsigned i = this->_size; i; i--)
    if (*contentsPtr++ == value)
      return TRUE;

  return FALSE;
}

template <class Type>
Boolean
SimpleArray<Type>::contains(Type value, unsigned start, unsigned end) const
{
  if ((end < start) || (end >= this->_size) || (start >= this->_size)) {
    std::cerr << "SimpleArray::contains called with invalid start (" << start
   << ") and end (" << end << ") arguments (array size: " << this->_size << ")" << std::endl;
    return FALSE;
  }

  register Type *contentsPtr = this->_contents + start;
  for (register unsigned i = end - start + 1; i; i--)
    if (*contentsPtr++ == value)
      return TRUE;

  return FALSE;
}

template <class Type>
Boolean
SimpleArray<Type>::containsOnly(Type value) const
{
  register Type *contentsPtr = this->_contents;
  for (register unsigned i = this->_size; i; i--)
    if (*contentsPtr++ != value)
      return FALSE;

  return TRUE;
}

template <class Type>
Boolean
SimpleArray<Type>::containsOnly(Type value, unsigned start, unsigned end) const
{
  if ((end < start) || (end >= this->_size) || (start >= this->_size)) {
    std::cerr << "SimpleArray::containsOnly called with invalid start (" << start
   << ") and end (" << end << ") arguments (array size: " << this->_size << ")" << std::endl;
    return FALSE;
  }

  register Type *contentsPtr = this->_contents + start;
  for (register unsigned i = end - start + 1; i; i--)
    if (*contentsPtr++ != value)
      return FALSE;

  return TRUE;
}

template <class Type> unsigned
SimpleArray<Type>::occurrencesOf(Type value, unsigned start, unsigned end) const
{
  if (end > this->_size - 1) {
    std::cerr << "Warning! SimpleArray::occurrencesOf() called with end=" << end
   << " on array of size " << this->_size << ". Truncated!" << std::endl;
    end = this->_size - 1;
  }

  if (start > end) {
    std::cerr << "Warning! SimpleArray::occurrencesOf() called with start > end"
   << std::endl;
    return 0;
  }

  unsigned N = 0;
  this->resetIterator(start);

  for (unsigned i = end - start + 1; i; i--)
    if ((*this)++ == value)
      N++;

  return N;
}

template <class Type>
int
SimpleArray<Type>::indexOf(Type value, int dir, unsigned start) const
{
  this->resetIterator(start);

  if (dir > 0) {
    for (register int i = this->_size - start; i; i--)
      if ((*this)++ == value)
  return this->_itIndex - 1;
  }
  else {
    for (register int i = start + 1; i; i--)
      if ((*this)-- == value)
  return this->_itIndex + 1;
  }

  return -1;
}

template <class Type>
void
SimpleArray<Type>::removeAll(Type value)
{
  if (!this->_size)
    return;

  unsigned i,j;
  for (i = 0, j = 0; i < this->_size; i++) {
    Type el = this->getElConst(i);
    if (el != value) {
      if (i != j)
        this->setEl(j, el);
      j++;
    }
  }

  this->newSize(j);
}

template <class Type>
void
SimpleArray<Type>::removeAllIn(Type floor, Type ceil, unsigned *N)
{
  if (!this->_size)
    return;

  if (floor == ceil)
    removeAll(floor);

  if (floor > ceil)
    std::swap(floor, ceil);

  unsigned i, j;
  unsigned n = 0;
  for (i = 0, j = 0; i < this->_size; i++) {
    Type value = this->getElConst(i);
    if ((value >= floor) || (value <= ceil))
      n++;
    else {
      if (i != j)
        this->setEl(j, value);
      j++;
    }
  }

  this->newSize(j);

  if (N)
    *N = n;
}

template <class Type>
void
SimpleArray<Type>::removeAllNotIn(Type floor, Type ceil,
          unsigned *nBelow, unsigned *nAbove)
{
  if (!this->_size)
    return;

  if (floor > ceil)
    std::swap(floor, ceil);

  unsigned belowCtr = 0;
  unsigned aboveCtr = 0;

  unsigned i,j;
  for (i = 0, j = 0; i < this->_size; i++) {
    Type value = this->getElConst(i);
    if (value < floor)
      belowCtr++;
    else if (value > ceil)
      aboveCtr++;
    else {
      if (i != j)
        this->setEl(j, value);
      j++;
    }
  }

  this->newSize(j);

  if (nAbove)
    *nAbove = aboveCtr;

  if (nBelow)
    *nBelow = belowCtr;
}

template <class Type>
UnsignedArray
SimpleArray<Type>::indicesOf(Type value) const
{
  UnsignedArray indices(0);

  this->resetIterator();
  for (unsigned i = 0; i < this->_size; i++)
    if ((*this)++ == value)
      indices.append(i);

  return indices;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::common(const SimpleArray<Type>& array) const
{
  SimpleArray<Type> common;

  Type *elementPtr = this->_contents;
  for (unsigned i = this->_size; i; i--) {
    if (array.contains(*elementPtr) && (!common.contains(*elementPtr)))
      common.append(*elementPtr);
    elementPtr++;
  }

  return common;
}

// #if HAVE_FINITE
// #ifndef finite
// extern "C" int finite(double);
// #endif /* finite() not defined (as macro) */
#define FINITE(x) finite(x)
// #elif HAVE_ISFINITE
// #ifndef isfinite
// extern "C" int isfinite(double);
// #endif /* isfinite() not defined (as macro) */
// #define FINITE(x) isfinite(x)
// #else
// #error "Neither finite() nor isfinite() is defined on your system"
// #endif /* HAVE_ISFINITE */


template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::prune()
{
  unsigned i,j;
  for (i = 0, j = 0; i < this->_size; i++) {
    double value = double(this->getElConst(i));
    if (FINITE(value)) {
      if (i != j)
        this->setEl(j, Type(value));
      j++;
    }
  }

  this->newSize(j);

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::randuniform(double min, double max)
{
  double range = max - min;

  for (unsigned i = 0; i < this->_size; i++)
    this->setEl(i, Type(drand48() * range + min));

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::randnormal(double mean, double std)
{
  for (unsigned i = 0; i < this->_size; i++)
    this->setEl(i, Type(gauss(mean, std)));

  return *this;
}

template <class Type>
UnsignedArray
SimpleArray<Type>::qsortIndexAscending() const
{
  unsigned i;

  if (!this->_size)
    return UnsignedArray(0);

  Type *elementPtr = this->_contents;
  IndexStruct<Type> *indexStructArray = new IndexStruct<Type>[this->_size];
  IndexStruct<Type> *indexStructPtr   = indexStructArray;
  for (i = 0; i < this->_size; i++) {
    indexStructPtr->element = *elementPtr++;
    (indexStructPtr++)->index = i;
  }

  ::qsort(indexStructArray, this->_size, sizeof(IndexStruct<Type>),
    IndexStruct<Type>::compareAscending);

  UnsignedArray sortedIndices(this->_size);
  unsigned     *sortedIndicesPtr = sortedIndices.contents();
  indexStructPtr = indexStructArray;
  for (i = 0; i < this->_size; i++)
    *sortedIndicesPtr++ = (indexStructPtr++)->index;

  delete [] indexStructArray;

  return sortedIndices;
}

template <class Type>
UnsignedArray
SimpleArray<Type>::qsortIndexDescending() const
{
  unsigned i;

  if (!this->_size)
    return UnsignedArray(0);

  Type *elementPtr = this->_contents;
  IndexStruct<Type> *indexStructArray = new IndexStruct<Type>[this->_size];
  IndexStruct<Type> *indexStructPtr   = indexStructArray;
  for (i = 0; i < this->_size; i++) {
    indexStructPtr->element = *elementPtr++;
    (indexStructPtr++)->index   = i;
  }

  ::qsort(indexStructArray, this->_size, sizeof(IndexStruct<Type>),
    IndexStruct<Type>::compareDescending);

  UnsignedArray sortedIndices(this->_size);
  unsigned     *sortedIndicesPtr = sortedIndices.contents();
  indexStructPtr = indexStructArray;
  for (i = 0; i < this->_size; i++)
    *sortedIndicesPtr++ = (indexStructPtr++)->index;

  delete [] indexStructArray;

  return sortedIndices;
}

template <class Type>
Type
SimpleArray<Type>::min(unsigned *index) const
{
  assert(this->_size);

  this->resetIterator();
  Type min = (*this)++;
  if (index)
    *index = 0;

  for (unsigned i = 1; i < this->_size; i++) {
    Type value = (*this)++;
    if (value < min) {
      min = value;
      if (index)
  *index = i;
    }
  }

  return min;
}

template <class Type>
Type
SimpleArray<Type>::max(unsigned *index) const
{
  assert(this->_size);

  this->resetIterator();
  Type max = (*this)++;
  if (index)
    *index = 0;

  for (unsigned i = 1; i < this->_size; i++) {
    Type value = (*this)++;
    if (value > max) {
      max = value;
      if (index)
  *index = i;
    }
  }

  return max;
}

template <class Type>
void
SimpleArray<Type>::extrema(Type *min, Type *max) const
{
  assert(this->_size);

  this->resetIterator();

  *max = *min = (*this)++;

//   if (this->_debug)
//       std::cout << this->_size << " :: " << *max << " :: " << *min << std::endl;

  for (unsigned i = 1; i < this->_size; i++) {
    Type value = (*this)++;
    if (value < *min)
      *min = value;
    if (value > *max)
      *max = value;
  }

//   if (this->_debug)
//       std::cout << this->_size << " :: " << *max << " :: " << *min << std::endl;
}

template <class Type>
Type
SimpleArray<Type>::range(unsigned *minIndex, unsigned *maxIndex) const
{
  assert(this->_size);

  this->resetIterator();
  Type min = (*this)++;
  if (minIndex)
    *minIndex = 0;

  Type max = min;
  if (maxIndex)
    *maxIndex = 0;

  for (unsigned i = 1; i < this->_size; i++) {
    Type value = (*this)++;
    if (value < min) {
      min = value;
      if (minIndex)
  *minIndex = i;
    }

    if (value > max) {
      max = value;
      if (maxIndex)
  *maxIndex = i;
    }
  }

  return (max - min);
}

template <class Type>
double
SimpleArray<Type>::sum() const
{
  double sum = 0;

  this->resetIterator();
  for (unsigned i = this->_size; i; i--)
    sum += asDouble((*this)++);

  return sum;
}

template <class Type>
double
SimpleArray<Type>::sum2() const
{
  double sum2 = 0;

  this->resetIterator();
  for (unsigned i = this->_size; i; i--) {
    double value = asDouble((*this)++);
    sum2 += value*value;
  }

  return sum2;
}

template <class Type>
double
SimpleArray<Type>::prod() const
{
  double prod = 0;

  if (this->_size) {
    this->resetIterator();
    prod = asDouble((*this)++);
    for (unsigned i = this->_size - 1; i; i--)
      prod *= asDouble((*this)++);
  }

  return prod;
}

template <class Type>
double
SimpleArray<Type>::prod2() const
{
  double prod2 = 0;

  if (this->_size) {
    this->resetIterator();
    prod2 = asDouble((*this)++);
    prod2 *= prod2;
    for (unsigned i = this->_size - 1; i; i--) {
      double value = asDouble((*this)++);
      prod2 *= value*value;
    }
  }

  return prod2;
}

template <class Type>
double
SimpleArray<Type>::var() const
{
  if (!this->_size)
    return 0;

  double sum    = 0;
  double sum2 = 0;

  this->resetIterator();
  for (unsigned i = this->_size; i; i--) {
    double value = asDouble((*this)++);
    sum    += value;
    sum2 += value*value;
  }

  return sum2/this->_size - SQR(sum/this->_size);
}

template <class Type>
Type
SimpleArray<Type>::median() const
{
  assert(this->_size);

  SimpleArray<Type> array(*this);

  return array.medianVolatile();
}

template <class Type>
Type
SimpleArray<Type>::medianVolatile()
{
  assert(this->_size);

  return _randomizedSelect(0, this->_size - 1,
                           (this->_size % 2) ? (this->_size + 1) / 2 : this->_size / 2);
}

template <class Type>
Type
SimpleArray<Type>::mode(const Type) const
{
  std::cerr << "Warning! SimpleArray::mode called but not implemented; returning median" << std::endl;

  return median();
}

template <class Type>
DblArray
SimpleArray<Type>::cumSum() const
{
  DblArray result(this->_size);

  if (!this->_size)
    return result;

  this->resetIterator();
  result.resetIterator();

  double value = asDouble((*this)++);
  result++ = value;
  for (unsigned i = this->_size - 1; i; i--) {
    value += asDouble((*this)++);
    result++ = value;
  }

  return result;
}

template <class Type>
DblArray
SimpleArray<Type>::cumProd() const
{
  DblArray result(this->_size);

  if (!this->_size)
    return result;

  double value = asDouble((*this)++);
  result++ = value;
  for (unsigned i = this->_size - 1; i; i--) {
    value *= asDouble((*this)++);
    result++ = value;
  }

  return result;
}

//
// Set functions
//

template <class Type>
void
SimpleArray<Type>::ceil(Type ceil)
{
  this->resetIterator();
  for (register unsigned i = 0; i < this->_size; i++)
    if ((*this)++ > ceil)
      this->setEl(i, ceil);
}

template <class Type>
void
SimpleArray<Type>::floor(Type floor)
{
  this->resetIterator();
  for (register unsigned i = 0; i < this->_size; i++)
    if ((*this)++ < floor)
      this->setEl(i, floor);
}

//
// Boolean operations
//

template <class Type>
Boolean
SimpleArray<Type>::operator != (const SimpleArray<Type>& array) const
{
  if (this->_size != array._size)
    return TRUE;

  this->resetIterator();
  array.resetIterator();

  for (unsigned i = this->_size; i; i--)
    if ((*this)++ != array++)
      return TRUE;

  return FALSE;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator && (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--, value1ptr++, value2ptr++)
      *boolPtr++ = *value1ptr && *value2ptr;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator || (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--, value1ptr++, value2ptr++)
      *boolPtr++ = *value1ptr || *value2ptr;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator == (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ == static_cast<Type>(value);
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator != (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ != static_cast<Type>(value);
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator >= (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ >= value;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator > (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ > value;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator <= (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ <= value;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator < (double value) const
{
  BoolArray boolArray(this->_size);

  if (this->_size) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *valuePtr = this->_contents;
    for (register unsigned i = this->_size; i; i--)
      *boolPtr++ = *valuePtr++ < value;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator >= (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--)
      *boolPtr++ = *value1ptr++ >= *value2ptr++;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator > (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--)
      *boolPtr++ = *value1ptr++ > *value2ptr++;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator <= (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--)
      *boolPtr++ = *value1ptr++ <= *value2ptr++;
  }

  return boolArray;
}

template <class Type>
BoolArray
SimpleArray<Type>::operator < (const SimpleArray<Type>& array) const
{
  BoolArray boolArray((Boolean) FALSE, this->_size);

  unsigned nElements = MIN(this->_size, array._size);
  if (nElements) {
    Boolean *boolPtr  = boolArray.contents();
    Type    *value1ptr = this->_contents;
    Type    *value2ptr = array._contents;
    for (register unsigned i = nElements; i; i--)
      *boolPtr++ = *value1ptr++ < *value2ptr++;
  }

  return boolArray;
}

//
// Arithmetic operations
//

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator += (Type value)
{
  this->resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ += value;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator += (const SimpleArray<Type>& array)
{
  assert(this->_size == array._size);

  this->resetIterator();
  array.resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ += array++;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator -= (Type value)
{
  this->resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ -= value;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator -= (const SimpleArray<Type>& array)
{
  assert(this->_size == array._size);

  this->resetIterator();
  array.resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ -= array++;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator *= (double value)
{
  this->resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ *= (Type)value;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator *= (const SimpleArray<Type>& array)
{
  assert(this->_size == array._size);

  this->resetIterator();
  array.resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ *= array++;

  return *this;
}

template <class Type>
SimpleArray<Type>&
SimpleArray<Type>::operator /= (const SimpleArray<Type>& array)
{
  assert(this->_size == array._size);

  this->resetIterator();
  array.resetIterator();
  for (unsigned i = this->_size; i; i--)
    (*this)++ /= array++;

  return *this;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::abs() const
{
  SimpleArray<Type> result(this->_size);

  Type   *sourcePtr = this->_contents;
  Type   *destPtr   = result.contents();
  for (unsigned i = this->_size; i; i--) {
    *destPtr++ = (*sourcePtr < Type(0)) ? Type(0)-(*sourcePtr) : *sourcePtr;
    sourcePtr++;
  }

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::round(unsigned n) const
{
  SimpleArray<Type> result(this->_size);

  Type   *sourcePtr = this->_contents;
  Type   *destPtr   = result.contents();
  if (n) {
    unsigned factor = (unsigned) pow(10.0, double(n));
    for (unsigned i = this->_size; i; i--) {
      *destPtr++ = ROUND(factor*double(*sourcePtr))/factor;
      sourcePtr++;
    }
  }
  else
    for (unsigned i = this->_size; i; i--) {
      *destPtr++ = ROUND(double(*sourcePtr));
      sourcePtr++;
    }

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::sqr() const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *resultPtr = result.contents();

  for (unsigned i = this->_size; i; i--, sourcePtr++)
    *resultPtr++ = SQR(*sourcePtr);

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::sqrt() const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *resultPtr = result.contents();

  for (unsigned i = this->_size; i; i--, sourcePtr++)
    *resultPtr++ = Type(::sqrt(asDouble(*sourcePtr))); // (bert) fix casting

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::operator ^ (int power) const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *resultPtr = result.contents();

  for (unsigned i = this->_size; i; i--)
    *resultPtr++ = Type(intPower(int(asDouble(*sourcePtr++)), power));

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::ln() const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *resultPtr = result.contents();

  for (unsigned i = this->_size; i; i--)
    *resultPtr++ = Type(::log(asDouble(*sourcePtr++))); // (bert) fix casting

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::log() const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *resultPtr = result.contents();

  for (unsigned i = this->_size; i; i--)
    *resultPtr++ = Type(::log10(asDouble(*sourcePtr++)));

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::exp() const
{
  return ::exp(1.0)^(*this);
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::exp10() const
{
  return 10^(*this);
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::sample(unsigned maxN) const
{
  double step = double(this->_size - 1)/(maxN - 1);

  if (step <= 1.0)
    return SimpleArray<Type>(*this);

  SimpleArray<Type> result(maxN);

  Type *destPtr = result._contents;
  double offset = 0;

  for (unsigned i = maxN; i; i--, offset += step)
    *destPtr++ = *(this->_contents + unsigned(::floor(offset)));

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::applyElementWise(Type (*function) (Type)) const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *destPtr   = result._contents;

  for (unsigned i = this->_size; i; i--)
    *destPtr++ = function(*sourcePtr++);

  return result;
}

template <class Type>
SimpleArray<Type>
SimpleArray<Type>::map(const ValueMap& map) const
{
  SimpleArray<Type> result(this->_size);

  Type *sourcePtr = this->_contents;
  Type *destPtr   = result._contents;

  for (unsigned i = this->_size; i; i--)
    *destPtr++ = Type(map(asDouble(*sourcePtr++)));

  return result;
}

//
// Type conversions
//

template <class Type>
SimpleArray<int>
asIntArray(const SimpleArray<Type>& array)
{
  IntArray cast(array.size());

  const Type *contentsPtr = array.contents();
  int  *castPtr     = cast.contents();

  for (register unsigned i = array.size(); i; i--)
    *castPtr++ = (int) *contentsPtr++;

  return cast;
}

template <class Type>
SimpleArray<float>
asFloatArray(const SimpleArray<Type>& array)
{
  SimpleArray<float> cast(array.size());

  const Type *contentsPtr = array.contents();
  float *castPtr = cast.contents();

  for (register unsigned i = array.size(); i; i--)
    *castPtr++ = (float) *contentsPtr++;

  return cast;
}

template <class Type>
SimpleArray<double>
asDblArray(const SimpleArray<Type>& array)
{
  DblArray cast(array.size());

  const Type *contentsPtr = array.contents();
  double     *castPtr     = cast.contents();

  for (register unsigned i = array.size(); i; i--)
    *castPtr++ = asDouble(*contentsPtr++);

  return cast;
}

//
// Friends
//

template <class Type>
SimpleArray<Type>
max(const SimpleArray<Type>& array1, const SimpleArray<Type>& array2)
{
  unsigned n = size(array1);

  if (size(array2) != n) {
    std::cerr << "max(SimpleArray, SimpleArray): arrays not of equal size" << std::endl;
    return SimpleArray<Type>(0);
  }

  SimpleArray<Type> result(array1);

  Type *resultPtr       = result.contents();
  const Type *array2Ptr = array2.contents();

  for (unsigned i = n; i; i--, resultPtr++, array2Ptr++)
    if (*array2Ptr > *resultPtr)
      *resultPtr = *array2Ptr;

  return result;
}

template <class Type>
SimpleArray<Type>
min(const SimpleArray<Type>& array1, const SimpleArray<Type>& array2)
{
  unsigned n = size(array1);

  if (size(array2) != n) {
    std::cerr << "min(SimpleArray, SimpleArray): arrays not of equal size" << std::endl;
    return SimpleArray<Type>(0);
  }

  SimpleArray<Type> result(array1);

  Type *resultPtr       = result.contents();
  const Type *array2Ptr = array2.contents();

  for (unsigned i = n; i; i--, resultPtr++, array2Ptr++)
    if (*array2Ptr < *resultPtr)
      *resultPtr = *array2Ptr;

  return result;
}

//
// Private functions
//

// Computes the i'th order statistic (the i'th smallest element) of
// an array.  This is the routine that does the real work, and is
// analogous to quicksort.

template <class Type>
Type
SimpleArray<Type>::_randomizedSelect(int p, int r, int i)
{
// If we're only looking at one element of the array, we must have
// found the i'th order statistic.

  if (p == r)
    return this->_contents[p];

// Partition the array slice A[p..r].  This rearranges the array so that
// all elements of A[p..q] are less than all elements of A[q+1..r].
// (Nothing fancy here, this is just a slight variation on the standard
// partition done by quicksort.)

  int q = _randomizedPartition(p, r);
  int k = q - p + 1;

// Now that we have partitioned A, its i'th order statistic is either
// the i'th order statistic of the lower partition A[p..q], or the
// (i-k)'th order stat. of the upper partition A[q+1..r] -- so
// recursively find it.

  if (i <= k)
    return _randomizedSelect(p, q, i);
  else
    return _randomizedSelect(q+1, r, i-k);
}

template <class Type>
int
SimpleArray<Type>::_randomizedPartition(int p, int r)
{
// Compute a random number between p and r

  //int i = (random() / (MAXINT / (r-p+1))) + p;
  // Changed because of problems locating random() on a Sun
  int i = int(ROUND(drand48() * (r-p+1) + p));

// Swap elements p and i

  std::swap(this->_contents[p], this->_contents[i]);

  return _partition(p, r);
}

template <class Type>
int
SimpleArray<Type>::_partition(int p, int r)
{
  Type x = this->_contents[p];
  int i  = p-1;
  int j  = r+1;

  while (1) {
    do { j--; } while (this->_contents[j] > x);
    do { i++; } while (this->_contents[i] < x);

    if (i < j)
      std::swap(this->_contents[i], this->_contents[j]);
    else
      return j;
  }
}

template <class Type>
SimpleArray<Type>
operator ^ (double base, const SimpleArray<Type>& array) {
  unsigned N = array.size();

  SimpleArray<Type> result(N);
  const Type *sourcePtr = array.contents();
  Type *resultPtr = result.contents();
  for (unsigned i = N; i != 0; i--)
    *resultPtr++ = Type(pow(base, double(*sourcePtr++)));

  return result;
}

// template <class Type> unsigned SimpleArray<Type>::_rangeErrorCount = 25;


#endif
