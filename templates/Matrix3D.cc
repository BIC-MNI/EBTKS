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
$RCSfile: Matrix3D.cc,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2003-04-16 15:09:37 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "EBTKS/Matrix3D.h"
#include "EBTKS/MatrixSupport.h"
#include "EBTKS/FileIO.h"
#include <sys/stat.h>

template <class Type> 
unsigned Mat3D<Type>::_rangeErrorCount = 10;

template <class Type> 
Boolean Mat3D<Type>::flushToDisk = TRUE;

//
// Some mathematical operations
//
template <class T1, class T2>
Mat3D<T1>& operator += (Mat3D<T1>& A, const Mat3D<T2>& B)
{
  unsigned nslis = A.getslis();

  if ((B.getslis() != nslis) || 
      (B.getrows() != A.getrows()) || (B.getcols() != A.getcols())) {
    cerr << "Matrices of incompatible sizes for +=" << endl;
    return A;
  }

  for(unsigned s = 0; s < nslis; s++)
    A[s] += B[s];
  
  return A;
}

template <class T1, class T2>
Mat3D<T1>& operator -= (Mat3D<T1>& A, const Mat3D<T2>& B)
{
  unsigned nslis = A.getslis();

  if ((B.getslis() != nslis) || 
      (B.getrows() != A.getrows()) || (B.getcols() != A.getcols())) {
    cerr << "Matrices of incompatible sizes for -=" << endl;
    return A;
  }

  for(unsigned s = 0; s < nslis; s++)
    A[s] -= B[s];
  
  return A;
}

template <class T1, class T2>
Mat3D<T1>& pmultEquals (Mat3D<T1>& A, const Mat3D<T2>& B)
{
  unsigned nslis = A.getslis();

  if ((B.getslis() != nslis) || 
      (B.getrows() != A.getrows()) || (B.getcols() != A.getcols())) {
    cerr << "Matrices of incompatible sizes for pmultEquals" << endl;
    return A;
  }

  for(unsigned s = 0; s < nslis; s++)
    pmultEquals(A[s], B[s]);
  
  return A;
}

template <class T1, class T2>
Mat3D<T1>& pdivEquals (Mat3D<T1>& A, const Mat3D<T2>& B)
{
  unsigned nslis = A.getslis();

  if ((B.getslis() != nslis) || 
      (B.getrows() != A.getrows()) || (B.getcols() != A.getcols())) {
    cerr << "Matrices of incompatible sizes for pdivEquals" << endl;
    return A;
  }

  for(unsigned s = 0; s < nslis; s++)
    pdivEquals(A[s], B[s]);
  
  return A;
}

/*
template <class T1, class T2>
Mat<T1> operator * (const Mat<T1>& A, const Mat<T2>& B)
{
  unsigned arows = A.getrows();
  unsigned acols = A.getcols();

  unsigned brows    = B.getrows();
  unsigned bcols    = B.getcols();
  unsigned bmaxcols = B.getmaxcols();

  Mat<T1> Temp(arows, bcols);

  if (acols != brows) {
    cerr << "Mat sizes incompatible for *" << endl;
    return Temp;
  }
  
  T1        *tempPtr = (T1 *) Temp.getEl()[0];
  const T1 **rowPtr  = A.getEl();
  
  for (unsigned i = arows; i; i--) {
    const T2 *argColPtr = B.getEl()[0];
    
    for (unsigned j = bcols; j; j--) {
      const T1 *ptr    = *rowPtr;
      const T2 *argPtr = argColPtr++;
      *tempPtr = (T1) 0;
      
      for (unsigned k = acols; k; k--) {
	*tempPtr += *ptr++ * (T1) *argPtr;
	argPtr += bmaxcols;
      }
      tempPtr++;
    }
    rowPtr++;
  }
  
  return Temp;
}
*/

/*********************************
Mat3D class definitions and member functions
**********************************/
//
//-------------------------// 
//
template <class Type>
Mat3D<Type>::Mat3D(const Mat3D<Type>& A)
{
  _slis = A._slis;  _rows = A._rows; _cols = A._cols;
  _slices = 0;
  _el = 0;
  
  _allocateEl();
  for(unsigned s=0; s<_slis ; s++)
    _slices[s] = A._slices[s];
  _setEl();
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>::Mat3D(unsigned nslis, unsigned nrows, unsigned ncols, Type value)
{
   if (((int)nslis < 0) || ((int)nrows < 0) || ((int)ncols < 0)) {
      cerr << "Error: cannot create Mat3D with negative dimension." << endl;
      exit(1);
   }

   _slis = nslis;
   _rows = nrows;
   _cols = ncols;
   _el = 0;
   _slices = 0;
   
   _allocateEl();
   
   if (value != 0){
      for(unsigned s=0; s<_slis; s++)
	 _slices[s].fill(value);
   }
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>::Mat3D(unsigned nslis, unsigned nrows, unsigned ncols)
{
   _slis = nslis;
   _rows = nrows;
   _cols = ncols;
   _el = 0;
   _slices = 0;
   
   _allocateEl();
}

//
//-------------------------// 
//
//this function will copy same block of data to every slice in image
template <class Type>
Mat3D<Type>::Mat3D(unsigned nslis, unsigned nrows, unsigned ncols,
			 const Type *data)
{
   _slis = nslis;
   _rows = nrows;
   _cols = ncols;
   _el = 0;
  _slices = 0; 

   _allocateEl();

   for(unsigned s=0; s<_slis; s++)
      _slices[s] = Mat<Type>(_rows, _cols, data);
   _setEl();
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>::Mat3D(const char *filename, int type)
{
  _slis = 0;
  _rows = 0;
  _cols = 0;
  _el = 0;
  _slices = 0;
  
  switch(type) {
  case RAW:
    loadRaw(filename);
    break;
  case ASCII:
    loadAscii(filename);
    break;
  }
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>::~Mat3D()
{
  clear();
}

////
template <class Type>
void
Mat3D<Type>::clear()
{
   if (_slices) {
      #ifdef DEBUG
      cout << "Freeing 3D Matrix " << _el << endl;
      #endif
      delete [] _slices;
      delete [] _el;
   }

   _slices = 0;
   _el = 0;
   _slis = _rows = _cols = 0;
}

//
//
//-------------------------// 
//

template <class Type>
Mat3D<Type>&
Mat3D<Type>::setSlice(unsigned slice, const Mat<Type>& A)
{
   if ((A.getrows() != _rows) || (A.getcols() != _cols)) {
      cerr << "Blabla" << endl;
      return *this;
   }

   if (slice >= _slis) {
      cerr << "Blabla" << endl;
      return *this;
   }

   _slices[slice] = A;
   _setEl();

   return *this;
}

/************************************************************************
Three dimensional Mat display functions
*************************************************************************/
  
/************************************************************************
 This display function displays each two dimensional slice of the entire
 matrix.  A seperate function is used for matrices with complex _elements
*************************************************************************/

template <class Type>
ostream&
Mat3D<Type>::display(ostream& os) const
{
  for(unsigned s=0; s < _slis; s++)
    os << s << endl << _slices[s] << endl;

  return os;
}

/************************************************************************
This display function displays the requested portion of the three
dimensional matrix.  The first two arguments are the start and stop
_slices(respectiv_ely)(first slice starts at zero) that are to be displayed.
The third and fourth arguments are the start and stop _rows and the fifth and
sixth arguments are the start and stop columns(all start at zero). A
seperate function is used for three dimensional matrices.
*************************************************************************/

template <class Type>
ostream&
Mat3D<Type>::display(ostream& os, unsigned s1, unsigned s2, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const
{
//Make sure that the stop dimensions are greater than the start dimensions
   
  if ((s1 > s2) || (r1 > r2) || (c1 > c2)){
    cerr << "Error in display: improper slice or row or column sizes." << endl;
    cerr << s1 << " to " << s2 << " and" << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }

//Make sure that the requested _slices are actually within the matrix
     
  if (s2 >= _slis) { 
    cerr<<"The requested _slices are not all defined for this matrix"<<endl;
    exit(1);
  }
    
  for(unsigned s=s1; s <= s2; s++)
    os << s << endl << _slices[s](r1,r2,c1,c2) << endl;

  return os;
}

/*****************************************************************************
Three dimensional matrix assignment operators******************************************************************************/

template <class Type>
Mat3D<Type>& 
Mat3D<Type>::operator = (const Mat3D<Type>& A)
{
  if (this == &A)
    return *this;

  if ((A._slis != _slis) || (A._rows != _rows) || (A._cols != _cols)) {
    _slis = A._slis;
    _rows = A._rows;
    _cols = A._cols;
    _allocateEl();
  }

  if (_slices && _el) {
    for(unsigned s=0; s < _slis ; s++)
      _slices[s] = A._slices[s];
    _setEl();
  }

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>& 
Mat3D<Type>::absorb(Mat3D<Type>& A)
{
  if (this == &A)
    return *this;

  // Delete current contents;
  if (_slices) {
    delete [] _slices;
    delete [] _el;
  }
  
  // Copy all from A
  _slis    = A._slis;
  _rows    = A._rows;
  _cols    = A._cols;
  _slices  = A._slices;
  _el      = A._el;
  
  // Empty A
  A._slis = A._rows = A._cols = 0;
  A._slices = 0;
  A._el = 0;

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Type
Mat3D<Type>::operator () (unsigned s, unsigned r, unsigned c) const 
{
  if ((s >= _slis) || (r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << s << ", " << r << ", " << c 
	   << ") exceed matrix dimensions. Changed to (" 
	   << ::min(s, _slis - 1) << ", " 
	   << ::min(r, _rows - 1) << ", " 
	   << ::min(c, _cols - 1) << ")" << endl;
      _rangeErrorCount--;
    }
    s = ::min(s, _slis - 1);
    r = ::min(r, _rows - 1);
    c = ::min(c, _cols - 1);
  }

  return(_el[s][r][c]);
}

//
//-------------------------// 
//
template <class Type>
Type& 
Mat3D<Type>::operator () (unsigned s, unsigned r, unsigned c) 
{
  if ((s >= _slis) || (r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << s << ", " << r << ", " << c 
	   << ") exceed matrix dimensions. Changed to (" 
	   << ::min(s, _slis - 1) << ", " 
	   << ::min(r, _rows - 1) << ", " 
	   << ::min(c, _cols - 1) << ")" << endl;
      _rangeErrorCount--;
    }
    s = ::min(s, _slis - 1);
    r = ::min(r, _rows - 1);
    c = ::min(c, _cols - 1);
  }

  return(_el[s][r][c]);
}

// 
//
// crop a section of the Mat3D
template <class Type>
Mat3D<Type>
Mat3D<Type>::operator () (unsigned s1, unsigned s2, unsigned r1, unsigned r2, 
			  unsigned c1, unsigned c2) const {
  if ((s1 > s2) || (s1 >= _slis) || (s2 >= _slis)) {
    cerr << "Error in cropting: improper _slis sizes." << endl;
    cerr << s1 << " to " << s2 << " and" << endl;
    exit(1);
  }
  
  Mat3D<Type> A(s2-s1+1,r2-r1+1,c2-c1+1);
  for(unsigned s=s1; s<= s2 ; s++)
    A.setSlice(s, _slices[s](r1,r2,c1,c2));
  return(A);     
}

/****************************************************************************
These overloaded operators allow accessing of individual _elements within the
3D matrix.  By convention, the first _element of the first slice of the 
matrix is lab_eled _element zero.  Since the matrix is stored in row major format, The _element numbers increase going to the right and at the following row. The _element at the last row and last column of the last slice corresponds to the largest possible value of the input argument.
*****************************************************************************/

/****************************************************************************
This accessing operator prevents accidental modification by having a constant
output argument. A seperate function is used for accessing _elements within
complex matrices.
*****************************************************************************/

template <class Type>
Type
Mat3D<Type>::operator () (unsigned n) const 
{
  if (n >= _slis*_rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _slis*_rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _slis*_rows*_cols - 1;
  }

  unsigned slice_index= n/(_rows*_cols);
  unsigned twoD_index= (n%(_rows*_cols));
  unsigned row_index= twoD_index/_cols;
  unsigned columnindex= twoD_index%_cols;   
  
  return _el[slice_index][row_index][columnindex]; 
}

/****************************************************************************
This accessing operator allows modification of the individual _elements of
the matrix.  A seperate function is used for accessing _elements within
complex matrices.
*****************************************************************************/
template <class Type>
Type& 
Mat3D<Type>::operator () (unsigned n) 
{
  if (n >= _slis*_rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _slis*_rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _slis*_rows*_cols - 1;
  }

  unsigned slice_index= n/(_rows*_cols);
  unsigned twoD_index= (n%(_rows*_cols));
  unsigned row_index= twoD_index/_cols;
  unsigned columnindex= twoD_index%_cols;   

  return _el[slice_index][row_index][columnindex];
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::operator () (const Mat3D<Type>& A)
{
   if ((A._slis != _slis) || (A._rows != _rows) || (A._cols != _cols)) {
      _slis = A._slis;
      _rows = A._rows;
      _cols = A._cols;
      _allocateEl();
   }
   
   for(unsigned s=0; s < _slis ; s++)
      _slices[s] = A._slices[s];
   _setEl();

   return *this;
}

//
//-------------------------// 
//
template <class Type>
int
Mat3D<Type>::operator != (const Mat3D<Type>& Arg) const
{
  if ((_slis != Arg._slis))
    return 1;
  
  for(unsigned s=0; s<_slis; s++)
    if (_slices[s] != Arg._slices[s])
      return 1;
  
  return 0;
}

template <class Type>
Mat3D<Type>& 
Mat3D<Type>::operator += (complex addend)
{
  for(unsigned s=0; s<_slis ; s++)
    _slices[s] += addend;

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::operator *= (complex scale)
{
   for(unsigned s=0; s<_slis ; s++)
      _slices[s] *= scale;
   return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::applyElementWise(double (*function)(double))
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].applyElementWise(function);
  
  return(*this);
}

template <class Type>
Mat3D<Type>&
Mat3D<Type>::applyElementWiseC2D(double (*function)(complex))
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].applyElementWiseC2D(function);
  
  return(*this);
}

template <class Type>
Mat3D<Type>&
Mat3D<Type>::applyElementWiseC2C(complex (*function)(complex))
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].applyElementWiseC2C(function);
  
  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::applyIndexFunction(IndexFunction3D F)
{
  for (unsigned s = 0; s < _slis; s++) {
    Type *elPtr = *(_el[s]);
    for (unsigned r = 0; r < _rows; r++)
      for (unsigned c = 0; c < _rows; c++)
	*elPtr++ = Type(F(s, r, c));
  }
  
  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::applyIndexFunction(ComplexIndexFunction3D F)
{
  for (unsigned s = 0; s < _slis; s++) {
    Type *elPtr = *(_el[s]);
    for (unsigned r = 0; r < _rows; r++)
      for (unsigned c = 0; c < _rows; c++)
	*elPtr++ = Type(F(s, r, c));
  }
  
  return(*this);
}

//
//-------------------------// 
//
template <class Type>
unsigned
Mat3D<Type>::length() const
{
  return ::max(_slis, _rows, _cols);
}

//
//-------------------------// 
//
template <class Type>
Type 
Mat3D<Type>::min(unsigned *sli, unsigned *row, unsigned *col) const
{
   Type min = ***_el;
   unsigned    sliOfMin = 0, rowOfMin = 0, colOfMin = 0;
   Type ***sliPtr = _el;
   for(unsigned s=_slis; s != 0; s--){
      Type **rowPtr = *sliPtr++;
      for(unsigned i = _rows ; i != 0 ; i--) {
	 Type *colPtr = *rowPtr++;
	 for(unsigned j = _cols ; j != 0 ; j--){
	    if (*colPtr < min) {
	       min = *colPtr;
	       rowOfMin = i;
	       colOfMin = j;
	       sliOfMin = s;
	    }
	    colPtr++;
	 }
      }
   }
   if (row)
      *row = rowOfMin;
   if (col)
      *col = colOfMin;
   if (sli)
      *sli = sliOfMin;
   
   return(min);
}

#ifdef USE_COMPMAT
complex
Mat3D<complex>::min(unsigned *, unsigned *, unsigned *) const
{
  cerr << "Mat3D<complex>::min() called but not implemented" << endl;
  return 0;
}
#endif

#ifdef USE_FCOMPMAT
fcomplex
Mat3D<fcomplex>::min(unsigned *, unsigned *, unsigned *) const
{
  cerr << "Mat3D<fcomplex>::min() called but not implemented" << endl;
  return 0;
}
#endif

//
//-------------------------// 
//
template <class Type>
Type 
Mat3D<Type>::max(unsigned *sli, unsigned *row, unsigned *col) const
{
   Type max = ***_el;
   unsigned   sliOfMax = 0, rowOfMax = 0, colOfMax = 0;
   
   for (unsigned s=0 ; s< _slis ; s++){
      Type **rowPtr = _el[s];
      for(unsigned i = _rows ; i != 0 ; i--) {
	 Type *colPtr = *rowPtr++;
	 for(unsigned j = _cols ; j != 0 ; j--){
	    if (*colPtr > max) {
	       max = *colPtr;
	       rowOfMax = i;
	       colOfMax = j;
	       colOfMax = s;
	    }
	    colPtr++;
	 }
      }
   }
   if (row)
      *row = rowOfMax;
   if (col)
      *col = colOfMax;
   if (sli)
      *sli = sliOfMax;
   
   return(max);
}

#ifdef USE_COMPMAT
complex
Mat3D<complex>::max(unsigned *, unsigned *, unsigned *) const
{
  cerr << "Mat3D<complex>::max() called but not implemented" << endl;
  return 0;
}
#endif

#ifdef USE_FCOMPMAT
fcomplex
Mat3D<fcomplex>::max(unsigned *, unsigned *, unsigned *) const
{
  cerr << "Mat3D<fcomplex>::max() called but not implemented" << endl;
  return 0;
}
#endif

//
//-------------------------// 
//
template <class Type>
Type
Mat3D<Type>::median(Type minVal, Type maxVal) const
{
  return asArray(minVal, maxVal).medianVolatile();
}

#ifdef USE_COMPMAT
complex
Mat3D<complex>::median(complex, complex) const
{
  cerr << "Mat3D<complex>::median() called but not implemented" << endl;
  return 0;
}
#endif

#ifdef USE_FCOMPMAT
fcomplex
Mat3D<fcomplex>::median(fcomplex, fcomplex) const
{
  cerr << "Mat3D<fcomplex>::median() called but not implemented" << endl;
  return 0;
}
#endif

//
//-------------------------// 
//
template <class Type>
complex
Mat3D<Type>::csum() const
{
  complex result = 0;
  for(unsigned s=0; s<_slis ; s++)
    result += _slices[s].csum();
  
  return result;
}

template <class Type>
complex
Mat3D<Type>::csum2() const
{
  complex result = 0;
  for(unsigned s=0; s<_slis ; s++)
    result += _slices[s].csum2();
  
  return result;
}

template <class Type>
complex
Mat3D<Type>::ctrace() const
{
   unsigned mindim = ::min(_slis, _rows, _cols);
   
   complex temp = 0;
   for (unsigned i=0 ; i < mindim ; i++)
      temp += _el[i][i][i];
   
   return(temp);
}

template <class Type>
complex
Mat3D<Type>::cdet() const
{
  cerr << "Mat3D::cdet() called but not implemented" << endl;
  return 0;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::exp()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].exp();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::log()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].log();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::cos()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].cos();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::sin()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].sin();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::abs()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].abs();
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::conj()
{
  for(unsigned s=0; s<_slis; s++)
    _slices[s].conj();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::round()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].round();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::sqrt()
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].sqrt();
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::pow(double exponent)
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].pow(exponent);
  return *this;
}

//
//-------------------------//
//
template < class Type>
Mat3D<Type>&
Mat3D<Type>::clip(Type minVal, Type maxVal, Type minFill, Type maxFill)
{
  for (unsigned s = 0; s < _slis; s++)
    _slices[s].clip(minVal, maxVal, minFill, maxFill);
  
  return(*this);
}

//
//-------------------------//
//
template < class Type>
Mat3D<Type>&
Mat3D<Type>::map(const ValueMap& valueMap)
{
  for(unsigned s = 0; s < _slis; s++)
    _slices[s].map(valueMap);

  return(*this);
}

//
//-------------------------// 
//
//The imagescale function will scale the image between 0 and 255
//unless given a range with a minimum and maximum value
//It also accept minimum and maximum values 
 
template < class Type>
Mat3D<Type>&
Mat3D<Type>::scale(double minout, double maxout, double minin, double maxin)
{
  if (minin >= maxin) {
    minin = asDouble(min());
    maxin = asDouble(max());
  }

  for(unsigned s = 0; s < _slis; s++)
    _slices[s].scale(minout, maxout, minin, maxin);
   
  return(*this);
}

//
//------------------------/
//
template <class Type>
Boolean
Mat3D<Type>::load(const char *filename, int type) const
{
  Boolean status = FALSE;

  switch(type) {
  case RAW:
    status = loadRaw(filename);
    break;
  case ASCII:
    status = loadAscii(filename);
    break;
  default:
    cerr << "Unrecognized type for loading" << endl;
    status = FALSE;
    break;
  }

  return status;
}

template <class Type>
Boolean
Mat3D<Type>::loadRaw(const char *filename, unsigned nslis, unsigned nrows, 
		     unsigned ncols)
{
  InputFile matrixFile(filename);
  
  if (!matrixFile) {
    cerr << "Error in loadRaw: error opening file." << endl;
    return FALSE;
  }
  
  _checkMatrixDimensions(filename, nslis, nrows, ncols);

  if ((nslis && (nslis != _slis)) || (nrows && (nrows != _rows)) || (_cols && (ncols != _cols))) {
    _slis = nslis;
    _rows = nrows;
    _cols = ncols;
    _allocateEl();
  }
  
  for(unsigned s=0; s<_slis; s++)
    if (!matrixFile.stream().read((unsigned char *)_el[s] + _rows*sizeof(Type *),
				  _rows*_cols*sizeof(Type)))
      return FALSE;

  return TRUE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat3D<Type>::loadAscii(const char *filename)
{
  InputFile infile(filename);
  
  if (!infile){
    cerr << "Error in loadAsccii: error opening file." << endl;
    return FALSE;
  }
  
  infile >> _slis >> _rows >> _cols;
  
  _allocateEl();
  
  for(unsigned s=0; s<_slis; s++) {
    for(unsigned i=0; i < _rows; i++) {
      for(unsigned j=0; j < _cols ; j++)
	infile >> _el[s][i][j];
    }
  }
  
  return TRUE;
}

#ifdef USE_COMPMAT
Boolean
Mat3D<complex>::loadAscii(const char *)
{
  cerr << "Mat3D<complex>::loadAscii() called but not implemented" << endl;
  return FALSE;
}
#endif

#ifdef USE_FCOMPMAT
Boolean
Mat3D<fcomplex>::loadAscii(const char *)
{
  cerr << "Mat3D<fcomplex>::loadAscii() called but not implemented" << endl;
  return FALSE;
}
#endif

template <class Type>
Boolean
Mat3D<Type>::save(const char *filename, int type) const
{
  Boolean status = FALSE;

  switch(type) {
  case RAW:
    status = saveRaw(filename);
    break;
  case ASCII:
    status = saveAscii(filename);
    break;
  default:
    cerr << "Unrecognized type for saving" << endl;
    status = FALSE;
    break;
  }

  return status;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat3D<Type>::saveRaw(const char *filename) const
{
  ofstream outfile(filename);
  if (!outfile) {
    cerr << "Error in saveRaw: error opening file." << endl;
    return FALSE;
  }

  Boolean status = TRUE;

  for(unsigned s=0; s<_slis && status; s++)
    if (!outfile.write((unsigned char *)_el[s] + _rows*sizeof(Type *), 
		       _rows*_cols*sizeof(Type)))
      status = FALSE;
    
  outfile.close();

  return status;
}

template <class Type>
Boolean
Mat3D<Type>::saveAscii(const char *filename) const
{
  ofstream outfile(filename);
  if (!outfile){
    cerr << "Error in saveAscciifile: error opening file." << endl;
    return FALSE;
  }

  Boolean status = TRUE;

  outfile << _slis << " " << _rows << " " << _cols << endl;

  if (!outfile)
    status = FALSE;
  
  for(unsigned s=0; s<_slis && status; s++) {
    for(unsigned i=0; i < _rows && status; i++) {
      for(unsigned j=0; j < _cols && status; j++)
	if (!(outfile << _el[s][i][j] << " "))
	  status = FALSE;
      outfile << endl;
    }
    outfile << endl;
  }
  
  outfile.close();
  
  return status;
}
//
//------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::insert(const Mat3D<Type>& A, int slice, int row, int col)
{
  int nSlis = A._slis;

  if ((slice + nSlis > _slis)) 
    cerr << "Warning: Mat3D<Type>::insert() in _slis" << endl;
    
  if (slice + nSlis > _slis)
      nSlis = _slis - slice;
    
  for(unsigned s=0; s< nSlis; s++)
      _slices[s+slice] = _slices[s+slice].insert(A._slices[s], row, col);
      
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::insert(const char *path, unsigned nslis, unsigned nrows, unsigned ncols, 
		  int slice, int row, int col)
{
  InputFile argFile(path);
  if (!argFile) {
    cerr << "Couldn't open file " << path << endl;
    return *this;
  }

  _checkMatrixDimensions(path, nslis, nrows, ncols);

  // Create a row-sized buffer
  Type *buffer = 0;
  allocateArray(ncols, buffer);
  if (!buffer) {
    cerr << "Couldn't allocate buffer" << endl;
    return *this;
  }

  int destSlice = slice;
  for (unsigned s = nslis; s; s--, destSlice++) {
    Boolean sliceValid = (destSlice >= 0) && (destSlice < _slis);
    int destRow    = row;
    for (unsigned i = nrows; i; i--, destRow++) {
      if (!(argFile.stream().read((unsigned char *) buffer, ncols*sizeof(Type)))) {
	cerr << "Error while reading file " << path << endl;
	freeArray(buffer);
	return *this;
      }
      if (sliceValid && (destRow >= 0) && (destRow < _rows)) {
	Type       *elPtr     = _el[destSlice][destRow] + col;
	const Type *bufferPtr = buffer;
	int         destCol   = col;
	
	for (unsigned j = ncols; j; j--, destCol++, elPtr++, bufferPtr++)
	  if ((destCol >= 0) && (destCol < _cols))
	    *elPtr = *bufferPtr;
      }
    }
  }

  freeArray(buffer);

  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::fill(Type value)
{
  for(unsigned s=0; s<_slis; s++)
    _slices[s].fill(value);
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::randuniform(double min, double max)
{
  for(unsigned s=0; s<_slis; s++)
    _slices[s].randuniform(min, max);
   
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::randnormal(double mean, double std)
{
  for(unsigned s=0; s<_slis; s++)
    _slices[s].randnormal(mean, std);
  
  return *this;
}
//
//-------------------------//
//
// Private functions start
//
template <class Type>
void
Mat3D<Type>::_allocateEl()
{
  if (_slices) {
#ifdef DEBUG
    cout << "Freeing 3D Matrix " << _el << endl;
#endif
    delete [] _slices;
    delete [] _el;
  }
  _slices = 0;
  _el = 0;

  if (_slis && _rows && _cols) {
    _slices = new Mat<Type> [_slis];
    assert(_slices);
  
    _el = new Type** [_slis];
    assert(_el);
  
#ifdef DEBUG
    cout << "Allocating 3D real matrix of " << _slis << " _slices at " << _el << endl;
#endif
    
    for (unsigned slice = 0; slice < _slis; slice++)
      _slices[slice].resize(_rows, _cols); //done this way to save memory
 
    _setEl();
  }
}

//
//----------------//
//
template <class Type>
void
Mat3D<Type>::_setEl()
{
  if (_slices && _el)
    for (unsigned s = 0; s < _slis; s++)
      _el[s] = (Type **) _slices[s].getEl();
}

template <class Type>
void
Mat3D<Type>::_checkMatrixDimensions(const char *path, unsigned& nslis,
				    unsigned& nrows, unsigned& ncols) const
{
  if (Path(path).hasCompressedExtension())
    return;
  
  struct stat buf;
  stat(path, &buf);
  inferDimensions(buf.st_size/sizeof(Type), nslis, nrows, ncols);
}

template <class Type>
Mat3D<Type>&
Mat3D<Type>::_fft(unsigned nslis, unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  // Verify dimensions of FFT
  if ((nslis > 1) && ((nslis < _slis) || !isPowerOf2(nslis))) {
    cerr << "Warning! Mat3D<Type>::fft():" << endl
	 << "  Requested # slices for FFT (" << nslis << ") invalid;" << endl;
    nslis = ::max(unsigned(4), nextPowerOf2(::max(nslis, _slis)));
    cerr << "  increased to " << nslis << endl;
  }

  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat3D<Type>::fft():" << endl
	 << "  Requested # rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = ::max(unsigned(4), nextPowerOf2(::max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat3D<Type>::fft():" << endl
	 << "  Requested # cols for FFT (" << ncols << ") invalid;" << endl;
    ncols = ::max(unsigned(4), nextPowerOf2(::max(ncols, _cols)));
    cerr << "  increased to " << ncols << endl;
  }

  Boolean doX = (ncols != 1) && (_cols != 1) ? TRUE : FALSE;
  Boolean doY = (nrows != 1) && (_rows != 1) ? TRUE : FALSE;
  Boolean doZ = (nslis != 1) && (_slis != 1) ? TRUE : FALSE;

  if (!nslis)
    nslis = ::max(unsigned(4), nextPowerOf2(_slis));
  else if (nslis == 1)
    nslis = _slis;

  if (!nrows)
    nrows = ::max(unsigned(4), nextPowerOf2(_rows));
  else if (nrows == 1)
    nrows = _rows;

  if (!ncols)
    ncols = ::max(unsigned(4), nextPowerOf2(_cols));
  else if (ncols == 1)
    ncols = _cols;

  // Pad matrix to the final FFT dimensions
  pad(nslis, nrows, ncols, (nslis - _slis)/2, (nrows - _rows)/2, (ncols - _cols)/2, 0);

  // Allocate temporary arrays
  double *real = 0;
  double *imag = 0;
  unsigned maxDim = ::max((doZ ? _slis : 1), (doY ? _rows : 1), (doX ? _cols : 1));

  allocateArray(maxDim, real);
  allocateArray(maxDim, imag);
  assert(real && imag);

  unsigned slice, row, col;

  // Take 1D FFT in X (row) direction
  if (doX)
    for (slice = 0; slice < _slis; slice++)
      for(row = 0; row < _rows; row++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	Type    *sourcePtr = _el[slice][row];
	for(col = _cols; col; col--) {
	  complex value(*sourcePtr++);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_cols, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el[slice][row];
	for(col = _cols; col; col--)
	  *sourcePtr++ = Type(abs(complex(*realPtr++, *imagPtr++)));
      }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for (slice = 0; slice < _slis; slice++)
      for(col = 0; col < _cols; col++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	Type    *sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  complex value(*sourcePtr);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	  sourcePtr += _cols;
	}
	
	// Calculate 1D FFT
	fftFunc(_rows, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  *sourcePtr = Type(abs(complex(*realPtr++, *imagPtr++)));
	  sourcePtr += _cols;
	}
      }

  // Take 1D FFT in Z (slice) direction
  if (doZ) {
    unsigned offset = 0;
    for (row = 0; row < _rows; row++)
      for (col = 0; col < _cols; col++, offset++) {
	// Fill temporary FFT array
	double    *realPtr   = real;
	double    *imagPtr   = imag;
	Type    ***sourcePtr = _el;
	for(slice = _slis; slice; slice--) {
	  complex value(*(**sourcePtr++ + offset));
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_slis, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el;
	for(slice = _slis; slice; slice--)
	  *(**sourcePtr++ + offset) = Type(abs(complex(*realPtr++, *imagPtr++)));
      }
  }
  
  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}

#ifdef USE_COMPMAT
Mat3D<complex>&
Mat3D<complex>::_fft(unsigned nslis, unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  // Verify dimensions of FFT
  if ((nslis > 1) && ((nslis < _slis) || !isPowerOf2(nslis))) {
    cerr << "Warning! Mat3D<complex>::fft():" << endl
	 << "  Requested # slices for FFT (" << nslis << ") invalid;" << endl;
    nslis = ::max(unsigned(4), nextPowerOf2(::max(nslis, _slis)));
    cerr << "  increased to " << nslis << endl;
  }

  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat3D<complex>::fft():" << endl
	 << "  Requested # rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = ::max(unsigned(4), nextPowerOf2(::max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat3D<complex>::fft():" << endl
	 << "  Requested # cols for FFT (" << ncols << ") invalid;" << endl;
    ncols = ::max(unsigned(4), nextPowerOf2(::max(ncols, _cols)));
    cerr << "  increased to " << ncols << endl;
  }

  Boolean doX = (ncols != 1) && (_cols != 1) ? TRUE : FALSE;
  Boolean doY = (nrows != 1) && (_rows != 1) ? TRUE : FALSE;
  Boolean doZ = (nslis != 1) && (_slis != 1) ? TRUE : FALSE;

  if (!nslis)
    nslis = ::max(unsigned(4), nextPowerOf2(_slis));
  else if (nslis == 1)
    nslis = _slis;

  if (!nrows)
    nrows = ::max(unsigned(4), nextPowerOf2(_rows));
  else if (nrows == 1)
    nrows = _rows;

  if (!ncols)
    ncols = ::max(unsigned(4), nextPowerOf2(_cols));
  else if (ncols == 1)
    ncols = _cols;

  // Pad matrix to the final FFT dimensions
  pad(nslis, nrows, ncols, (nslis - _slis)/2, (nrows - _rows)/2, (ncols - _cols)/2, 0);

  // Allocate temporary arrays
  double *real = 0;
  double *imag = 0;
  unsigned maxDim = ::max((doZ ? _slis : 1), (doY ? _rows : 1), (doX ? _cols : 1));

  allocateArray(maxDim, real);
  allocateArray(maxDim, imag);
  assert(real && imag);

  unsigned slice, row, col;

  // Take 1D FFT in X (row) direction
  if (doX)
    for (slice = 0; slice < _slis; slice++)
      for(row = 0; row < _rows; row++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	complex *sourcePtr = _el[slice][row];
	for(col = _cols; col; col--) {
	  complex value(*sourcePtr++);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_cols, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el[slice][row];
	for(col = _cols; col; col--)
	  *sourcePtr++ = complex(*realPtr++, *imagPtr++);
      }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for (slice = 0; slice < _slis; slice++)
      for(col = 0; col < _cols; col++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	complex *sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  complex value(*sourcePtr);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	  sourcePtr += _cols;
	}
	
	// Calculate 1D FFT
	fftFunc(_rows, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  *sourcePtr = complex(*realPtr++, *imagPtr++);
	  sourcePtr += _cols;
	}
      }

  // Take 1D FFT in Z (slice) direction
  if (doZ) {
    unsigned offset = 0;
    for (row = 0; row < _rows; row++)
      for (col = 0; col < _cols; col++, offset++) {
	// Fill temporary FFT array
	double    *realPtr   = real;
	double    *imagPtr   = imag;
	complex ***sourcePtr = _el;
	for(slice = _slis; slice; slice--) {
	  complex value(*(**sourcePtr++ + offset));
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_slis, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el;
	for(slice = _slis; slice; slice--)
	  *(**sourcePtr++ + offset) = complex(*realPtr++, *imagPtr++);
      }
  }
  
  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}
#endif

#ifdef USE_FCOMPMAT
Mat3D<fcomplex>&
Mat3D<fcomplex>::_fft(unsigned nslis, unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  // Verify dimensions of FFT
  if ((nslis > 1) && ((nslis < _slis) || !isPowerOf2(nslis))) {
    cerr << "Warning! Mat3D<fcomplex>::fft():" << endl
	 << "  Requested # slices for FFT (" << nslis << ") invalid;" << endl;
    nslis = ::max(unsigned(4), nextPowerOf2(::max(nslis, _slis)));
    cerr << "  increased to " << nslis << endl;
  }

  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat3D<fcomplex>::fft():" << endl
	 << "  Requested # rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = ::max(unsigned(4), nextPowerOf2(::max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat3D<fcomplex>::fft():" << endl
	 << "  Requested # cols for FFT (" << ncols << ") invalid;" << endl;
    ncols = ::max(unsigned(4), nextPowerOf2(::max(ncols, _cols)));
    cerr << "  increased to " << ncols << endl;
  }

  Boolean doX = (ncols != 1) && (_cols != 1) ? TRUE : FALSE;
  Boolean doY = (nrows != 1) && (_rows != 1) ? TRUE : FALSE;
  Boolean doZ = (nslis != 1) && (_slis != 1) ? TRUE : FALSE;

  if (!nslis)
    nslis = ::max(unsigned(4), nextPowerOf2(_slis));
  else if (nslis == 1)
    nslis = _slis;

  if (!nrows)
    nrows = ::max(unsigned(4), nextPowerOf2(_rows));
  else if (nrows == 1)
    nrows = _rows;

  if (!ncols)
    ncols = ::max(unsigned(4), nextPowerOf2(_cols));
  else if (ncols == 1)
    ncols = _cols;

  // Pad matrix to the final FFT dimensions
  pad(nslis, nrows, ncols, (nslis - _slis)/2, (nrows - _rows)/2, (ncols - _cols)/2, 0);

  // Allocate temporary arrays
  double *real = 0;
  double *imag = 0;
  unsigned maxDim = ::max((doZ ? _slis : 1), (doY ? _rows : 1), (doX ? _cols : 1));

  allocateArray(maxDim, real);
  allocateArray(maxDim, imag);
  assert(real && imag);

  unsigned slice, row, col;

  // Take 1D FFT in X (row) direction
  if (doX)
    for (slice = 0; slice < _slis; slice++)
      for(row = 0; row < _rows; row++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	fcomplex *sourcePtr = _el[slice][row];
	for(col = _cols; col; col--) {
	  fcomplex value(*sourcePtr++);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_cols, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el[slice][row];
	for(col = _cols; col; col--)
	  *sourcePtr++ = fcomplex(*realPtr++, *imagPtr++);
      }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for (slice = 0; slice < _slis; slice++)
      for(col = 0; col < _cols; col++) {
	// Fill temporary FFT array
	double  *realPtr   = real;
	double  *imagPtr   = imag;
	fcomplex *sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  fcomplex value(*sourcePtr);
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	  sourcePtr += _cols;
	}
	
	// Calculate 1D FFT
	fftFunc(_rows, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = *_el[slice] + col;
	for(row = _rows; row; row--) {
	  *sourcePtr = fcomplex(*realPtr++, *imagPtr++);
	  sourcePtr += _cols;
	}
      }

  // Take 1D FFT in Z (slice) direction
  if (doZ) {
    unsigned offset = 0;
    for (row = 0; row < _rows; row++)
      for (col = 0; col < _cols; col++, offset++) {
	// Fill temporary FFT array
	double    *realPtr   = real;
	double    *imagPtr   = imag;
	fcomplex ***sourcePtr = _el;
	for(slice = _slis; slice; slice--) {
	  fcomplex value(*(**sourcePtr++ + offset));
	  *realPtr++ = ::real(value);
	  *imagPtr++ = ::imag(value);
	}
	
	// Calculate 1D FFT
	fftFunc(_slis, real, imag);
	
	// Put results back
	realPtr   = real;
	imagPtr   = imag;
	sourcePtr = _el;
	for(slice = _slis; slice; slice--)
	  *(**sourcePtr++ + offset) = fcomplex(*realPtr++, *imagPtr++);
      }
  }
  
  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}
#endif

template <class Type>
ostream&
operator << (ostream& os, const Mat3D<Type>& A)
{
  return A.display(os);
}

//***********************pad function ******************************

template <class Type>
Mat3D<Type>&
Mat3D<Type>::pad(unsigned nslis, unsigned nrows, unsigned ncols,
		 int slice, int row, int col, Type value)
{
  if ((nslis == _slis) && (nrows == _rows) && (ncols == _cols) && 
      !slice && !row && !col)
    return *this;

  char tempFile[256];
  get_temp_filename(tempFile);

  // Save the current matrix to disk
  if (flushToDisk && saveRaw(tempFile)) {
    unsigned argSlis = _slis;
    unsigned argRows = _rows;
    unsigned argCols = _cols;

    clear();

    _slis = nslis;
    _rows = nrows;
    _cols = ncols;

    _allocateEl();
    fill(value);

    insert(tempFile, argSlis, argRows, argCols, slice, row, col);
  }
  else {
    // Saving fails, create a new matrix
    Mat3D<Type> result(nslis, nrows, ncols, value);
    result.insert(*this, slice, row, col);
    
    this->absorb(result);
  }

  unlink(tempFile);

  return *this;
}

//
//-----------------------//
//

/***************************morphology functions*******************************/
//note works for odd and even se assumung the even se
//is centered at the right lowest pix_el of the four center pix_els.
#ifdef USE_DBLMAT
template <class Type>
Mat3D<Type>
Mat3D<Type>::erode(const Mat3D<double>& strel) const
{
  unsigned strelSlices  = strel.getslis();
  unsigned strelHeight = strel.getrows();
  unsigned strelWidth  = strel.getcols();

  if (((strelHeight == 1) && (strelWidth == 1) && (strelSlices == 1)) || !strelHeight || !strelWidth || !strelSlices)
    return Mat3D<Type>(*this);

  Mat3D<Type> padMatrix(pad(strelSlices/2, strelHeight/2, strelWidth/2));
  Mat3D<Type> result(_slis, _rows, _cols);
  
  unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
  unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
  unsigned n;
  
  for (unsigned s = 0; s < _slis ; s++){
     int padk=0;
     Type *resultPtr = result._el[s][0];

     for (unsigned y = _rows; y != 0; y--) {
        for (unsigned x = _cols; x != 0; x--) {
           n=s;
           double  minimum = MAXDOUBLE;
        
           for (register unsigned k = 0; k < strelSlices; k++) {
              Type *padMatrixPtr2 = padMatrix._el[n++][0] + padk;
              const double *strelPtr = strel.getEl()[k][0];
              for (register unsigned j = strelHeight; j != 0; j--) {
	         for (register unsigned i = strelWidth; i != 0; i--) {
                    if (*strelPtr >= 0)
                    minimum = MIN(minimum, (double) *padMatrixPtr2 - *strelPtr);
                    strelPtr++;
	            padMatrixPtr2++;
	         }
	         padMatrixPtr2 += padMatrixPtr2incr;
              }
           }
           *resultPtr = Type(minimum);
           resultPtr++;
           padk++;
        }
          padk += padMatrixPtr1incr;
     }   

  }  //for (unsigned s ...
  
  return result;
}

#ifdef USE_COMPMAT
Mat3D<complex>
Mat3D<complex>::erode(const Mat3D<double>&) const
{
  cerr << "Mat3D<complex>::erode() called but not implemented" << endl;
  return Mat3D<complex>(*this);
}
#endif

#ifdef USE_FCOMPMAT
Mat3D<fcomplex>
Mat3D<fcomplex>::erode(const Mat3D<double>&) const
{
  cerr << "Mat3D<fcomplex>::erode() called but not implemented" << endl;
  return Mat3D<fcomplex>(*this);
}
#endif

template <class Type>
Mat3D<Type>
Mat3D<Type>::dilate(const Mat3D<double>& strel) const
{
  unsigned strelSlices  = strel.getslis();
  unsigned strelHeight = strel.getrows();
  unsigned strelWidth  = strel.getcols();
  unsigned sl=0;
  unsigned r=0;
  unsigned c=0;

  if (((strelHeight == 1) && (strelWidth == 1) && (strelSlices == 1)) || !strelHeight || !strelWidth || !strelSlices)
    return Mat3D<Type>(*this);

   if ((strelSlices % 2) == 0) {
      strelSlices++;
      sl++;
   }

   if ((strelHeight % 2) == 0) {
      strelHeight++;
      r++;
   }
   
   if ((strelWidth % 2) == 0) {
      strelWidth++;
      c++;
   }

//   strel.display(cout);
//   cout << endl;
   
   Mat3D<double> newstrel(strelSlices, strelHeight, strelWidth, -1);

//   newstrel.display(cout);
//   cout << endl;
   
   newstrel.insert(strel.rotate180(), sl, r, c);
   
//   newstrel.display(cout);
//   cout << endl;

  Mat3D<Type> padMatrix(pad(strelSlices/2, strelHeight/2, strelWidth/2));
  Mat3D<Type> result(_slis, _rows, _cols);

  unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
  unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
  unsigned n;
  
  for (unsigned s = 0; s < _slis ; s++){
     int padk=0;
     Type *resultPtr = result._el[s][0];
     
     for (unsigned y = _rows; y != 0; y--) {
        for (unsigned x = _cols; x != 0; x--) {
           n=s;
           double  maximum = -MAXDOUBLE;
	   
           for (register unsigned k = 0; k < strelSlices; k++) {
              Type *padMatrixPtr2 = padMatrix._el[n++][0] + padk;
              const double *strelPtr = newstrel.getEl()[k][0];
              for (register unsigned j = strelHeight; j != 0; j--) {
	         for (register unsigned i = strelWidth; i != 0; i--) {
		    if (*strelPtr >= 0)
		       maximum = MAX(maximum, (double) *padMatrixPtr2 + *strelPtr);
		    strelPtr++;
	            padMatrixPtr2++;
	         }
	         padMatrixPtr2 += padMatrixPtr2incr;
              }
           }
           *resultPtr = Type(maximum);
           resultPtr++;
           padk++;
        }
          padk += padMatrixPtr1incr;
     }   

  }  //for (unsigned s ...
  
  return result;
}
#ifdef USE_COMPMAT
Mat3D<complex>
Mat3D<complex>::dilate(const Mat3D<double>&) const
{
  cerr << "Mat3D<complex>::dilate() called but not implemented" << endl;
  return Mat3D<complex>(*this);
}
#endif

#ifdef USE_FCOMPMAT
Mat3D<fcomplex>
Mat3D<fcomplex>::dilate(const Mat3D<double>&) const
{
  cerr << "Mat3D<fcomplex>::dilate() called but not implemented" << endl;
  return Mat3D<fcomplex>(*this);
}
#endif

#endif

//************************
//reverse the filter or rotate by 180
//used for convolution and morphology 
//-------------------------// 
//
template <class Type>
Mat3D<Type> 
Mat3D<Type>::rotate180() const           
{
   Mat3D<Type> Temp(_slis,_rows,_cols);
   unsigned ks,kr,kc;
   for (unsigned s=0 ; s < _slis ; s++) { 
      ks = _slis -s -1;  
      for (unsigned i=0 ; i < _rows ; i++) {
         kr = _rows -i -1;
         for (unsigned j=0 ; j < _cols ; j++){
	    kc = _cols -j -1;
	    Temp(ks,kr,kc) = _el[s][i][j];
         }
      }
   }
   
   return Temp;
}

//3D histogram
//this histogram works in the following way
//if there is no input in the histogram ex: histogram(), then
//it takes the defaults 0,0 and does the histogram from the smallest
//value in the volume to the largest value.
//if the input is histogram(0,255) it will give the histogram
//between 0 and 255.  If the input is 10,20 it will give only 10,20
//values


//Can be written in 3D return or 2D return

//Mat3D<Type>

/*template <class Type>
Mat<Type>
Mat3D<Type>::histogram(int  minin, int  maxin) const
{
   if(minin == maxin){
      minin=(int)min();
      maxin=(int)max();
    }
   int range = maxin - minin + 1;
   Mat<Type> histogram1(1,range);
   //Mat3D<Type> histogram1(1,1,range);
   Type *histPtr = (Type *)histogram1.getEl()[0];
   //Type *histPtr = histogram1._el[0][0];

   Type ***sliPtr = _el;
   for (unsigned s = _slis; s != 0; s--) {
      Type **rowPtr = *sliPtr++;
      for (unsigned i = _rows; i != 0; i--) {
	 Type *colPtr = *rowPtr++;
	 for (unsigned j = _cols; j != 0; j--) {
	    if ((*colPtr >= minin) && (*colPtr <= maxin)) {
	       int roundedEl = (int) rint((double)*colPtr);
	       (*(histPtr + roundedEl - minin))++;
	       //(*(histPtr + roundedEl - minin))++;
	    }
	    colPtr++;
	 }
      }

   }

return(histogram1);

}
*/

template <class Type>
Histogram
Mat3D<Type>::histogram(double minin, double maxin, unsigned n) const
{
  if (maxin <= minin) {
    minin = min();
    maxin = max();
  }

  Histogram histogram(minin, maxin, n);

  for(unsigned s=0; s < _slis; s++)
    histogram += _slices[s].histogram(minin, maxin, n);
  
  return(histogram);
}

#ifdef USE_COMPMAT
Histogram
Mat3D<complex>::histogram(double, double, unsigned) const
{
  cerr << "Mat3D<complex>::histogram() called but not implemented" << endl;
  return Histogram(0);
}
#endif

#ifdef USE_FCOMPMAT
Histogram
Mat3D<fcomplex>::histogram(double, double, unsigned) const
{
  cerr << "Mat3D<fcomplex>::histogram() called but not implemented" << endl;
  return Histogram(0);
}
#endif

//
//-------------------------//
//
template <class Type>
Mat3D<Type>&
Mat3D<Type>::histmod(const Histogram& hist1, const Histogram& hist2)
{
  return map(equalize(hist1, hist2));
}

//
//-------------------------//
//
template <class Type>
SimpleArray<Type>
Mat3D<Type>::asArray(Type minVal, Type maxVal) const
{
  unsigned N = 0;
  Type  ***slisPtr;
  Type   **rowPtr;
  Type    *elPtr;

  if (maxVal <= minVal) {
    minVal = min();
    maxVal = max();
    N = nElements();
  }
  else {
    // Scan matrix to determine # elements in the range
    slisPtr = _el;
    for (unsigned s = _slis; s; s--) {
      rowPtr = *slisPtr++;
      for (unsigned i = _rows; i; i--) {
	elPtr = *rowPtr++;
	for (unsigned j = _cols; j; j--, elPtr++)
	  if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	    N++;
      }
    }
  }
  
  SimpleArray<Type> array(N);

  if (N) {
    Type *arrayPtr = array.contents();
    slisPtr = _el;
    for (unsigned s = _slis; s; s--) {
      rowPtr = *slisPtr++;
      for (unsigned i = _rows; i; i--) {
	elPtr = *rowPtr++;
	for (unsigned j = _cols; j; j--)
	  if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	    *arrayPtr++ = *elPtr++;
      }
    }
  }
    
  return array;
}

template <class Type>
Mat3D<Type>&
Mat3D<Type>::ifft(unsigned nslis, unsigned nrows, unsigned ncols)
{
  _fft(nslis, nrows, ncols, ::ifft);
  unsigned factor = 1;
  if (nslis != 1) factor *= _slis;
  if (nrows != 1) factor *= _rows;
  if (ncols != 1) factor *= _cols;
  
  return *this /= factor;
}

// 
// Type conversions
//

#ifdef USE_DBLMAT
template <class Type> DblMat3D asDblMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  DblMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asDblMat(A[s]));

  return cast;
}
#endif

///
#ifdef USE_FLMAT
template <class Type> FlMat3D asFlMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  FlMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asFlMat(A[s]));

  return cast;
}
#endif

///
#ifdef USE_INTMAT
template <class Type> IntMat3D asIntMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  IntMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asIntMat(A[s]));

  return cast;
}
#endif

///
#ifdef USE_UINTMAT
template <class Type> UIntMat3D asUIntMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  UIntMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asUIntMat(A[s]));

  return cast;
}
#endif

///
#ifdef USE_SHMAT
template <class Type> ShMat3D asShMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  ShMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asShMat(A[s]));
  return cast;
}
#endif

///
#ifdef USE_USHMAT
template <class Type> UShMat3D asUShMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  UShMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asUShMat(A[s]));
  return cast;
}
#endif

///
#ifdef USE_CHRMAT
template <class Type> ChrMat3D asChrMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  ChrMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asChrMat(A[s]));
  return cast;
}
#endif

///
#ifdef USE_UCHRMAT
template <class Type> UChrMat3D asUChrMat(const Mat3D<Type>& A)
{
  unsigned _slis = A.getslis();
  UChrMat3D  cast(_slis, A.getrows(), A.getcols());
  for(unsigned s=0; s < _slis; s++)
    cast.setSlice(s, asUChrMat(A[s]));
  return cast;
}
#endif

//
//-------------------------// 
//
#ifdef USE_COMPMAT
template <class Type> 
Mat3D<complex> asCompMat(const Mat3D<Type>& A)
{
  unsigned nslis = A.getslis();

  Mat3D<complex> cast(nslis, A.getrows(), A.getcols());

  for (unsigned s = 0; s < nslis; s++)
    cast.setSlice(s, asCompMat(A[s]));

  return cast;
}

template <class Type> 
Mat3D<complex> asCompMat(const Mat3D<Type>& Re, const Mat3D<Type>& Im)
{
  unsigned nslis = Re.getslis();

  if ((Im.getslis() != nslis) || 
      (Im.getcols() != Re.getcols()) || (Im.getrows() != Re.getrows())) {
    cerr << "asCompMat: Re and Im matrices don't have the same dimensions; using Re only"
	 << endl;
    return asCompMat(Re);
  }

  Mat3D<complex> cast(nslis, Re.getrows(), Re.getcols());
  
  for (unsigned s = 0; s < nslis; s++)
    cast.setSlice(s, asCompMat(Re[s], Im[s]));

  return cast;
}
#endif

#ifdef USE_FCOMPMAT
template <class Type> 
Mat3D<fcomplex> asFcompMat(const Mat3D<Type>& A)
{
  unsigned nslis = A.getslis();
  Mat3D<fcomplex> cast(nslis, A.getrows(), A.getcols());

  for (unsigned s = 0; s < nslis; s++)
    cast.setSlice(s, asFcompMat(A[s]));

  return cast;
}

template <class Type> 
Mat3D<fcomplex> asFcompMat(const Mat3D<Type>& Re, const Mat3D<Type>& Im)
{
  unsigned nslis = Re.getslis();
  if ((Im.getslis() != nslis) || 
      (Im.getcols() != Re.getcols()) || (Im.getrows() != Re.getrows())) {
    cerr << "asFcompMat: Re and Im matrices don't have the same dimensions; using Re only"
	 << endl;
    return asFcompMat(Re);
  }

  Mat3D<fcomplex> cast(nslis, Re.getrows(), Re.getcols());
  
  for (unsigned s = 0; s < nslis; s++)
    cast.setSlice(s, asFcompMat(Re[s], Im[s]));

  return cast;
}
#endif

#ifdef USE_COMPMAT
#ifdef USE_DBLMAT
Mat3D<double>
applyElementWiseC2D(const Mat3D<complex>& A, double (*function)(const complex&))
{
  Mat3D<double> T(A._slis, A._rows, A._cols);

  for (unsigned s = 0; s < A._slis; s++)
    T[s] = applyElementWiseC2D(A[s], function);

  return T;
}

Mat3D<double>
arg(const Mat3D<complex>& A)
{
  unsigned nslis = A.getslis();

  Mat3D<double> T(nslis, A.getrows(), A.getcols());

  for (unsigned s = 0; s < nslis; s++)
    T.setSlice(s, arg(A[s]));

  return T;
}

Mat3D<double>
real(const Mat3D<complex>& A)
{
  return applyElementWiseC2D(A, real);
}

Mat3D<double>
imag(const Mat3D<complex>& A)
{
  return applyElementWiseC2D(A, imag);
}
#endif
#endif

#ifdef USE_FCOMPMAT
#ifdef USE_FLMAT
Mat3D<float>
applyElementWiseC2D(const Mat3D<fcomplex>& A, double (*function)(const complex&))
{
  Mat3D<float> T(A._slis, A._rows, A._cols);

  for (unsigned s = 0; s < A._slis; s++)
    T[s] = applyElementWiseC2D(A[s], function);

  return T;
}

Mat3D<float>
arg(const Mat3D<fcomplex>& A)
{
  unsigned nslis = A.getslis();

  Mat3D<float> T(nslis, A.getrows(), A.getcols());

  for (unsigned s = 0; s < nslis; s++)
    T.setSlice(s, arg(A[s]));

  return T;
}

Mat3D<float>
real(const Mat3D<fcomplex>& A)
{
  return applyElementWiseC2D(A, real);
}

Mat3D<float>
imag(const Mat3D<fcomplex>& A)
{
  return applyElementWiseC2D(A, imag);
}
#endif
#endif

//////////////////////////////////////////////////////////////////////
