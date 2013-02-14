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
$RCSfile: Array.h,v $
$Revision: 1.4 $
$Author: stever $
$Date: 2003-11-17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
/*
 * Essentially the template Array class described in Lippman's C++ Primer
 *
 * 11/12/1993
 *
 * Alex Zijdenbos
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <assert.h>
#include <string.h>

#include "EBTKS/trivials.h"
#include "EBTKS/MTypes.h"


const int DEFAULT_SIZE   = 0;
const int SIZE_INCREMENT = 32;

const int BEFORE = -1;
const int AFTER  = 1;

/******************
 * Array base class
 ******************/

template <class Type>
class Array {
// data members
  Array<Type> *_self; // Provided to circumvent the const mechanism for the iterator

protected:
  static unsigned int _arrayCtr;
  static Boolean  _debug;
  static unsigned int _rangeErrorCount;

  unsigned  _size;
  unsigned  _maxSize;
  Type     *_contents;

  // Iterator member
  unsigned  _itIndex;

public:
  static unsigned debug(Boolean on = TRUE);

// Manager functions
  Array (unsigned sz = DEFAULT_SIZE);
  Array (const Type&, unsigned sz);
  Array (const Type *, unsigned);
  Array (const Array&);
  virtual ~Array ();

  // Access functions
  inline Type&       operator [] (unsigned i);
  inline Type&       operator [] (int i);
  inline const Type& operator [] (unsigned i) const;
  inline const Type& operator [] (int i) const;

  inline virtual Type& getEl(unsigned i);
  inline virtual const Type& getEl(unsigned i) const;
  inline virtual const Type& getElConst(unsigned i) const;
  inline virtual void  setEl(unsigned i, Type value);

  // Iterator functions
  //virtual void  resetIterator(unsigned i = 0);       // Resets iterator to element i
  virtual void  resetIterator(unsigned i = 0) const;
  inline virtual Type& current();
  inline virtual const Type& current() const;

  // Prefix ascending iterators
  inline virtual Type& operator ++();
  inline virtual const Type&  operator ++() const;
  // Postfix ascending iterators
  inline virtual Type& operator ++(int);
  inline virtual const Type& operator ++(int) const;

  // Prefix descending iterators
  inline virtual Type& operator --();
  inline virtual const Type& operator --() const;
  // Postfix descending iterators
  inline virtual Type& operator --(int);
  inline virtual const Type& operator --(int) const;

  Array  operator () (unsigned nElements) const;           // Sub-array (first nElements)
  Array  operator () (unsigned start, unsigned end) const; // Sub-array

  Boolean operator ! () const;
  
  virtual Type&    first();
  virtual Type&    last();
  virtual unsigned size() const;
  // Direct access to contents - volatile!!
  virtual const Type *contents() const;
  virtual Type *contents();
  virtual operator const Type *() const;
  virtual operator Type *();

// Other functions
  Array&  operator = (const Array&);        // Copy
  Array&  absorb(Array&);                   // Copy (source destroyed)
  Array&  operator () (const Type *, unsigned); // Copy from C array
  // Returns C array. Uses <array> if specified; otherwise allocates a new array.
  Type   *asCarray(Type *array = 0) const;  
  Array&  append(const Type value);        // Append a value at the end
  Array&  append(const Array&);             // Concatenate two arrays
  Array&  insert(const Type&, unsigned index = 0);  // Insert an element at <index>
  Array&  insert(const Array&, unsigned index = 0); // Insert an array at <index>
  Array&  replace(const Array&, unsigned index = 0);// Replace, starting at <index>
  Type    remove(unsigned index = 0);       // Remove the element at <index>
  Type    removeLast();                     // Remove the last element
  Array&  reorder(const Array<unsigned>& indices); // Reorder the array
  Array&  shuffle();                        // Randomly re-order the array
  virtual void clear(const Type& value);
  virtual void newSize(unsigned);    // Changes size of array; keeps contents if possible
  Array&  destroy();                 // Destroys array (sets size to zero, frees memory)

  // Rotate operators
  Array& operator << (unsigned);
  Array& operator >> (unsigned);

  Array   sample(unsigned maxN) const;
  Array   applyElementWise(Type (*function) (Type)) const;

// I/O
  virtual std::ostream& print(std::ostream& os) const;

protected:
  void _grow(unsigned amount);       // Allocate more space; do not change _size
  virtual void _rangeError(unsigned& index) const;
  virtual void _notImplementedError() const;
};

template <class Type>
unsigned size(const Array<Type>& array);

template <class Type> 
std::ostream& operator << (std::ostream& os, const Array<Type>& array);

template <class Type>
Type&
Array<Type>::operator [] (unsigned i) 
{ 
  return getEl(i); 
}

template <class Type>
Type&
Array<Type>::operator [] (int i) 
{ 
  return getEl((unsigned)i); 
}

template <class Type>
const Type& 
Array<Type>::operator [] (unsigned i) const 
{ 
  return getElConst(i); 
}

template <class Type>
const Type& 
Array<Type>::operator [] (int i) const 
{ 
  return getElConst(i); 
}

//
// Access functions
//

template <class Type>
Type& 
Array<Type>::getEl(unsigned i) 
{
  if (i >= _size)_rangeError(i); 
  return _contents[i]; 
}

template <class Type>
const Type& 
Array<Type>::getEl(unsigned i) const
{
  return getElConst(i); 
}

template <class Type>
const Type& 
Array<Type>::getElConst(unsigned i) const
{ 
  if (i >= _size) _rangeError(i); 
  return _contents[i]; 
}
 
template <class Type>
void
Array<Type>::setEl(unsigned i, Type value)  
{ 
  if (i >= _size) _rangeError(i); 
  _contents[i] = value; 
}

template <class Type>
Type& 
Array<Type>::current()
{ 
  return _contents[_itIndex]; 
}

template <class Type>
const Type& 
Array<Type>::current() const  
{ 
  return _contents[_itIndex]; 
}

// Prefix ascending iterators
template <class Type>
Type& 
Array<Type>::operator ++()          
{ 
  return _contents[++_itIndex]; 
}

template <class Type>
const Type&  
Array<Type>::operator ++() const 
{
  return _contents[++_self->_itIndex]; 
}

// Postfix ascending iterators
template <class Type>
Type& 
Array<Type>::operator ++(int)       
{ 
  return _contents[_itIndex++]; 
}

template <class Type>
const Type& 
Array<Type>::operator ++(int) const 
{ 
  return _contents[_self->_itIndex++]; 
}

// Prefix descending iterators
template <class Type>
Type& 
Array<Type>::operator --()
{ 
  return _contents[--_itIndex]; 
}

template <class Type>
const Type& 
Array<Type>::operator --() const 
{ 
  return _contents[--_self->_itIndex]; 
}

// Postfix descending iterators
template <class Type>
Type& 
Array<Type>::operator --(int)
{ 
  return _contents[_itIndex--]; 
}

template <class Type>
const Type& 
Array<Type>::operator --(int) const 
{ 
  return _contents[_self->_itIndex--]; 
}


template <class Type>
unsigned
Array<Type>::debug(Boolean on) 
{
  unsigned N = _arrayCtr;

  if (_debug = on)
    _arrayCtr = 0;
    
  return N;
}

//
// constructors for Array class
//

template <class Type>
Array<Type>::Array (unsigned sz)
{
  _self = this;
  _size = _maxSize = sz;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    std::cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << std::flush;
  }
}

template <class Type>
Array<Type>::Array (const Type& value, unsigned sz)
{
  _self = this;
  _size = _maxSize = sz;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
    clear(value);
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    std::cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << std::flush;
  }
}

template <class Type>
Array<Type>::Array(const Type *initArray, unsigned nElements)
{
  _self = this;
  _size = _maxSize = nElements;

  if (_size) {
    _contents = new Type[_size];
    assert(_contents);
    memcpy(_contents, initArray, _size*sizeof(Type));
  }
  else
    _contents = 0;

  if (_debug) {
    _arrayCtr++;
    std::cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << std::flush;
  }
}

//
// Copy constructor
//

template <class Type>
Array<Type>::Array (const Array<Type>& array)
{
  _self     = this;
  _size     = _maxSize = 0;
  _contents = 0;

  operator = (array);

  if (_debug) {
    _arrayCtr++;
    std::cout << "C" << _arrayCtr << ":" << long(this) << ":" << _size << " " << std::flush;
  }
}

//
// destructor
//

template <class Type>
Array<Type>::~Array ()
{
  if (_debug) {
    _arrayCtr--;
    std::cout << "D" << _arrayCtr << ":" << long(this) << ":" << _size << " " << std::flush;
  }

  destroy();
}

//
// Iterator functions
//
/*
template <class Type>
void
Array<Type>::resetIterator(unsigned i)
{
  if (!_size)
    return;

  assert(i < _size);
  _itIndex = i;

  return;
}
*/
template <class Type>
void
Array<Type>::resetIterator(unsigned i) const
{
  if (!_size)
    return;

  assert(i < _size);
  _self->_itIndex = i;

  return;
}

template <class Type> 
Array<Type>
Array<Type>::operator () (unsigned nElements) const
{
  if (nElements > _size) {
    if (_rangeErrorCount) {
      std::cerr << "Warning! Array::operator(" << nElements 
     << ") called with on array of size " << _size << ". Value truncated!" << std::endl;
      _rangeErrorCount--;
    }
    nElements = _size;
  }
  
  Array<Type> subArray(nElements);

  Type *sourcePtr = _contents;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = nElements; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type> Array<Type>
Array<Type>::operator () (unsigned start, unsigned end) const
{
  unsigned n = end - start + 1;
  if (start + n > _size) {
    if (_rangeErrorCount) {
      std::cerr << "Warning! Array::operator(" << start << ", " << end
     << ") called with on array of size " << _size << ". Truncated!" << std::endl;
      _rangeErrorCount--;
    }
    n = _size - start;
  }

  Array<Type> subArray(n);

  Type *sourcePtr = _contents + start;
  Type *destPtr   = subArray._contents;
  for (register unsigned i = n; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(subArray);
}

template <class Type>
Boolean 
Array<Type>::operator ! () const 
{
  return Boolean(!_size); 
}

template <class Type>
Type&
Array<Type>::first() 
{ 
  return _contents[0]; 
}

template <class Type>
Type&
Array<Type>::last()
{ 
  return _contents[_size - 1]; 
}

template <class Type> 
unsigned 
Array<Type>::size() const
{ 
  return _size; 
}

template <class Type> 
const Type *
Array<Type>::contents() const 
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Type *
Array<Type>::contents()
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>::operator const Type *() const 
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>::operator Type *()
{ 
  return (_size) ? _contents : 0; 
}

template <class Type> 
Array<Type>& 
Array<Type>::operator = (const Array<Type>& array ) 
{
  if (this == &array) return *this;

  newSize(array.size());

  resetIterator();
  array.resetIterator();

  for (unsigned i = _size; i; i--)
    (*this)++ = array++;

  return *this;
}

template <class Type> Array<Type>& 
Array<Type>::absorb(Array<Type>& array) {
  if (this == &array) return *this;

  if (_contents)
    delete [] _contents;
  _size = _maxSize = array._size;
  _contents = array._contents;

  array._size = 0;
  array._contents = 0;

  return(*this);
}

template <class Type> Array<Type>& 
Array<Type>::operator () (const Type *newContents, unsigned size) 
{
  if (size > _maxSize) {
    if (_contents)
      delete [] _contents;
    _contents = new Type[_size = _maxSize = size];
    assert(_contents);
  }
  else
    _size = size;

  register const Type *sourcePtr = newContents;
  register Type *destPtr   = _contents;
  for (register unsigned i = _size; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return(*this);
}

template <class Type> 
Type *
Array<Type>::asCarray(Type *destPtr) const
{
  if (!_size)
    return 0;

  if (!destPtr)
    destPtr = new Type[_size];

  if (destPtr) {
    const Type *sourcePtr = _contents;
    for (register unsigned i = _size; i; i--)
      *destPtr++ = *sourcePtr++;
  }

  return destPtr;
}

template <class Type>
Array<Type>&  
Array<Type>::append(const Type value) 
{
  if (_maxSize<=_size)_grow(SIZE_INCREMENT); 
  _contents[_size++]=value; 
  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::append(const Array<Type>& array)
{
  unsigned nToAdd = array._size;

  if (!nToAdd)
    return *this;

  unsigned oldSize = _size;

  newSize(_size + nToAdd);

  const Type *sourcePtr = array._contents;
  Type *destPtr         = _contents + oldSize;
  for (register unsigned i = nToAdd; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::insert(const Type& value, unsigned index)
{
  if (index > _size) {
    if (_rangeErrorCount) {
      std::cerr << "Warning! Attempt to insert element outside range of array" << std::endl;
      _rangeErrorCount--;
    }
    return *this;
  }
  
  if (index == _size)
    return append(value);
  
  if (_maxSize <= _size)
    _grow(SIZE_INCREMENT);

  register Type *sourcePtr = _contents + _size - 1;
  register Type *destPtr   = sourcePtr + 1;
  for (register unsigned i = _size - index; i != 0; i--)
    *destPtr-- = *sourcePtr--;
  *destPtr = value;

  _size++;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::insert(const Array<Type>& array, unsigned index)
{
  if (array._size == 0)
    return *this;

  unsigned oldSize = _size;

  newSize(_size + array._size);

  unsigned i;
  Type *sourcePtr = _contents + oldSize - 1;
  Type *destPtr   = sourcePtr + array._size;
  for (i = oldSize - index; i != 0; i--)
    *destPtr-- = *sourcePtr--;

  sourcePtr = array._contents + array._size - 1;
  for (i = array._size; i != 0; i--)
    *destPtr-- = *sourcePtr--;

  return *this;
}

template <class Type> 
Array<Type>& 
Array<Type>::replace(const Array<Type>& array, unsigned index)
{
  if (array._size == 0)
    return *this;

  if (index + array._size > _size)
    newSize(index + array._size);

  register Type *sourcePtr = array._contents;
  register Type *destPtr   = _contents + index;
  for (register unsigned i = array._size; i != 0; i--)
    *destPtr++ = *sourcePtr++;

  return *this;
}

template <class Type>
Type
Array<Type>::remove(unsigned index)
{
  if (!_size) {
    if (_rangeErrorCount) {
      _rangeErrorCount--;
      std::cerr << "Warning! Attempt to remove element from empty array" << std::endl;
    }
      
    return _contents[0];
  }

  if (index >= _size)
    _rangeError(index);

  if (index == _size - 1) {
    _size--;
    return _contents[index];
  }

  Type value(_contents[index]);

  register Type *destPtr   = _contents + index;
  register Type *sourcePtr = destPtr + 1;
  for (unsigned i = _size - index - 1; i != 0; i--)
    *destPtr++ = *sourcePtr++;
  _size--;

  return value;
}

template <class Type>
Type
Array<Type>::removeLast()
{
  if (!_size) {
    if (_rangeErrorCount) {
      _rangeErrorCount--;
      std::cerr << "Warning! Attempt to remove element from empty array" << std::endl;
    }
    return _contents[0];
  }
    
  _size--;
  
  return _contents[_size];
}

template <class Type>
Array<Type>&
Array<Type>::reorder(const Array<unsigned>& indices)
{
  Array<Type>     temp(*this);
  Type           *elPtr  = _contents;
  const unsigned *indexPtr = indices.contents();

  unsigned n = std::min(indices.size(), _size);

  for (unsigned i = n; i; i--, elPtr++, indexPtr++)
    if (*indexPtr < _size)
      *elPtr = temp[*indexPtr];

  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::shuffle()
{
  for (unsigned i = 0; i < _size; i++) {
    unsigned j = unsigned(drand48()*_size);
    if (i != j) {
      Type temp    = _contents[i];
      _contents[i] = _contents[j];
      _contents[j] = temp;
    }
  }

  return *this;
}

template <class Type>
void
Array<Type>::clear(const Type& value)
{
  resetIterator();
  for (unsigned i = _size; i; i--)
    (*this)++ = value;
}

template <class Type>
void
Array<Type>::newSize(unsigned size)
{
  if (size == _size)
    return;

  if (size <= _maxSize) {
    _size = size;
    return;
  }

  Type *newContents = new Type[size];
  assert(newContents);

  if (_size) {
    register Type *sourcePtr = _contents;
    register Type *destPtr   = newContents;
    for (register unsigned i = _size; i != 0; i--)
      *destPtr++ = *sourcePtr++;
  }

  if (_contents)
    delete [] _contents;

  _contents = newContents;
  _size = _maxSize = size;
}

template <class Type>
Array<Type>&
Array<Type>::destroy()
{
  if (_contents) {
    delete [] _contents; 
    _contents = 0;
  }

  _size = _maxSize = 0;

  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::operator >> (unsigned n)
{
  if (_size) {
    n %= _size;
    Array<Type> temp(n);

    unsigned i;
    Type *sourcePtr = _contents + _size - 1;
    Type *destPtr   = temp._contents + n - 1;
    for (i = n; i; i--)
      *destPtr-- = *sourcePtr--;

    destPtr = _contents + _size - 1;
    for (i = _size - n; i; i--)
      *destPtr-- = *sourcePtr--;

    sourcePtr = temp._contents + n - 1;
    for (i = n; i; i--)
      *destPtr-- = *sourcePtr--;
  }
  return *this;
}

template <class Type>
Array<Type>&
Array<Type>::operator << (unsigned n)
{
  if (_size) {
    n %= _size;
    Array<Type> temp(n);

    unsigned i;
    Type *sourcePtr = _contents;
    Type *destPtr   = temp._contents;
    for (i = n; i; i--)
      *destPtr++ = *sourcePtr++;

    destPtr = _contents;
    for (i = _size - n; i; i--)
      *destPtr++ = *sourcePtr++;

    sourcePtr = temp._contents;
    for (i = n; i; i--)
      *destPtr++ = *sourcePtr++;
  }
  return *this;
}

template <class Type>
Array<Type>
Array<Type>::sample(unsigned maxN) const
{
  double step = double(_size - 1)/(maxN - 1);

  if (step <= 1.0)
    return Array<Type>(*this);

  Array<Type> result(maxN);
  
  Type *destPtr = result._contents;
  double offset = 0;
    
  for (unsigned i = maxN; i; i--, offset += step)
    *destPtr++ = *(_contents + unsigned(floor(offset)));

  return result;
}

template <class Type>
Array<Type>
Array<Type>::applyElementWise(Type (*function) (Type)) const
{
  Array<Type> result(_size);

  Type *sourcePtr = _contents;
  Type *destPtr   = result._contents;

  for (unsigned i = _size; i != 0; i--)
    *destPtr++ = function(*sourcePtr++);

  return result;
}

template <class Type> 
std::ostream&
Array<Type>::print(std::ostream& os) const
{
  std::cerr << "Array<Type>::print(): Cannot print an Array" << std::endl;
  return os;
}

//
// Private functions
//

template <class Type> 
void
Array<Type>::_grow(unsigned amount)
{
  unsigned size = _size;
  newSize(_maxSize + amount);
  _size = size;
}

template <class Type> 
void
Array<Type>::_rangeError(unsigned& index) const
{
  if (_rangeErrorCount) {
    _rangeErrorCount--;
    std::cerr << "Corrected: index " << index << " into array of size " << _size 
   << " !" << std::endl;
  }

  index = size() - 1;
}

template <class Type> 
void
Array<Type>::_notImplementedError() const
{
  std::cerr << "Array function called but not implemented" << std::endl;
  assert(0);
}

template <class Type> 
unsigned 
size(const Array<Type>& array) 
{
  return array.size(); 
}

template <class Type> 
std::ostream& 
operator << (std::ostream& os, const Array<Type>& array) 
{ 
  return array.print(os); 
}

template <class Type> unsigned int Array<Type>::_arrayCtr = 0;
template <class Type> Boolean  Array<Type>::_debug = FALSE;
template <class Type> unsigned int Array<Type>::_rangeErrorCount = 25;


#endif
