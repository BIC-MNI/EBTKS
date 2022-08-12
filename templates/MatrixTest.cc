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
$RCSfile: MatrixTest.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:26 $
$State: Exp $
--------------------------------------------------------------------------*/
#include <config.h>
#include "EBTKS/Matrix.h"
#include "EBTKS/MatrixSupport.h"
//#include "FileIO.h"  //included to load compressed data
/*********************************************************
template <class Type> const Mat<Type>::MATLAB = 0;
template <class Type> const Mat<Type>::RAW = 1;
template <class Type> const Mat<Type>::ASCII  = 2;
//template <class Type>  tells belong to template class 
*********************************************************/

//Template
/*********************************
Mat class definitions and member functions
**********************************/
//
//-------------------------// 
//

template <class Type>
Mat<Type>::Mat(const Mat<Type>& A)
{
   rows = A.rows; _cols = A._cols;
   maxrows = A.maxrows; _maxcols = A._maxcols;
   _el = 0;
   
   allocateEl();
   
   memcpy(*_el, *A._el, maxrows*_maxcols*sizeof(Type));
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols, Type)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::~Mat()
{
   clear();
}

//
//-------------------------// 
//

template <class Type>
void
Mat<Type>::clear()
{
   if (_el) {
      #ifdef DEBUG
      cout << "Freeing allocated memory at " << _el << endl;
      #endif
      delete ((char *) _el); //free((char *) _el);
      _el = 0;
   }

   _rows = _maxrows = _cols = _maxcols = 0;
}

template <class Type>
void
Mat<Type>::_allocateEl()
{
   if (_el) {
#ifdef DEBUG
      cout << "Freeing allocated memory at " << _el << endl;
#endif
     delete( (char *) _el);
     //  free((char *) _el);
   }

   unsigned int nBytes = _maxrows*sizeof(Type *) + _maxrows*_maxcols*sizeof(Type);
 // _el = (Type **) calloc(nBytes, 1); 
  _el = (Type **) new char [nBytes];
   
   if (!_el) {
      cerr << endl << "Error making pointer for " << _maxrows << " x " << _maxcols 
	 << " matrix" << endl;
      exit(1);
   }
   memset(_el, 0, nBytes);  //set all elements to zero

#ifdef DEBUG
   cout << "Allocated " << nBytes << " bytes at " <<  _el << endl;
#endif
   
   _el[0] = (Type *)(&(_el[_maxrows]));
   for (unsigned i=1 ; i < _maxrows ; i++)
      _el[i] = (Type *)(&(_el[i-1][_maxcols]));
}

/*******************************************/
/*******************************************/
/*******************************************/
/*******************************************/
/*******************************************/
/*******************************************/

Mat<complex>::Mat(const Mat<complex>& A)
{
  _rows = A._rows; _cols = A._cols;
  _maxrows = A._maxrows; _maxcols = A._maxcols;
  _el = 0;
   
  _allocateEl();
   
  memcpy(*_el, *A._el, _maxrows*_maxcols*sizeof(complex));
}

//
//-------------------------// 
//

Mat<complex>::Mat(unsigned nrows, unsigned ncols, complex)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;

   _el = 0;
   
   _allocateEl();
}

//
//-------------------------// 
//

Mat<complex>::~Mat()
{
   clear();
}

//
//-------------------------// 
//

void
Mat<complex>::clear()
{
   if (_el) {
      #ifdef DEBUG
      cout << "Freeing allocated memory at " << _el << endl;
      #endif
      delete ((char *) _el); //free((char *) _el);
      _el = 0;
   }

   _rows = _maxrows = _cols = _maxcols = 0;
}

void
Mat<complex>::_allocateEl()
{
   if (_el) {
#ifdef DEBUG
      cout << "Freeing allocated memory at " << _el << endl;
#endif
     delete( (char *) _el);
     //  free((char *) _el);
   }

   unsigned int nBytes = _maxrows*sizeof(complex *) + _maxrows*_maxcols*sizeof(complex);
 // _el = (complex **) calloc(nBytes, 1); 
  _el = (complex **) new char [nBytes];
   
   if (!_el) {
      cerr << endl << "Error making pointer for " << _maxrows << " x " << _maxcols 
	 << " matrix" << endl;
      exit(1);
   }
   memset(_el, 0, nBytes);  //set all elements to zero

#ifdef DEBUG
   cout << "Allocated " << nBytes << " bytes at " <<  _el << endl;
#endif
   
   _el[0] = (complex *)(&(_el[_maxrows]));
   for (unsigned i=1 ; i < _maxrows ; i++)
      _el[i] = (complex *)(&(_el[i-1][_maxcols]));
}

