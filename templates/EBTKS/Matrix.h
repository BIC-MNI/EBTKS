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
$RCSfile: Matrix.h,v $
$Revision: 1.7 $
$Author: rotor $
$Date: 2010-05-18 23:01:21 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MATRIX_H
#define _MATRIX_H

#include "EBTKS/SimpleArray.h"
#include "EBTKS/FileIO.h"  //included to load compressed data

/******************************************************************************
 * The user must use the following #defines for each matrix type that is 
 * required. Selective use of these #defines will greatly reduce mushrooming
 * of object code due to (unneccesary) template instantiations. Although most
 * compilers provide facilities to control template instantiations, this
 * approach was found more robust, for a large part due to the type conversions
 * provided by this libary.
 *   Without using the #defines, basic functionality will still be available,
 * but the typedefs and some functionality (such as type conversions) will 
 * be missing.
 *
 * Alex Zijdenbos, Aug 11, 1995.
 *
 * #define     type              typedef provided
 *-----------------------------------------------
 * USE_UCHRMAT  (unsigned char)   UChrMat
 * USE_CHRMAT   (char)            ChrMat
 * USE_USHMAT   (unsigned short)  UShMat
 * USE_SHMAT    (short)           ShMat
 * USE_UINTMAT  (unsigned int)    UIntMat
 * USE_INTMAT   (int)             IntMat
 * USE_FLMAT    (float)           FlMat
 * USE_DBLMAT   (double)          DblMat
 * USE_COMPMAT  (dcomplex)        CompMat
 * USE_FCOMPMAT (fcomplex)        fCompMat
 *****************************************************************************/

#include "EBTKS/dcomplex.h"

#ifdef USE_FCOMPMAT
  #include "EBTKS/fcomplex.h"
#endif

#ifdef HAVE_MATLAB
extern "C" {
  #include"mat.h" 
}
#endif

/* #include <stream.h> (bert) removed */
#include <math.h>
#include <fstream>		/* (bert) changed from fstream.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>
#include <iostream>		/* (bert) changed from iostream.h */
#include "EBTKS/MTypes.h"
#include "EBTKS/MatrixSupport.h"
#include "EBTKS/Histogram.h"

//VF: fixing templates
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef MIN
#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))
#endif


//VF: fixing stupid habit of putting templates into .cc file....
using namespace std;

const int MATLAB = 0;
const int RAW = 1;
const int ASCII  = 2;

typedef double  (*IndexFunction)(unsigned r, unsigned c);
typedef dcomplex (*ComplexIndexFunction)(unsigned r, unsigned c);

/******************** Explicit type conversions *******************/

// Forward declarations
template <class Type> class Mat;
template <class Type> class MatrixIterator;
template <class Type> class ConstMatrixIterator;

#ifdef USE_DBLMAT
//Converts Mat<Type> to Mat<double>
typedef Mat<double> DblMat;
template <class Type> Mat<double> asDblMat(const Mat<Type>&);
#endif

#ifdef USE_FLMAT
//Converts Mat<Type> to Mat<float>
typedef Mat<float> FlMat;
template <class Type> Mat<float> asFlMat(const Mat<Type>&);
#endif

#ifdef USE_INTMAT
//Converts Mat<Type> to Mat<int>
typedef Mat<int> IntMat;
template <class Type> Mat<int> asIntMat(const Mat<Type>&);
#endif

#ifdef USE_UINTMAT
//Converts Mat<Type> to Mat<unsigned int>
typedef Mat<unsigned int> UIntMat;
template <class Type> Mat<unsigned int> asUIntMat(const Mat<Type>&);
#endif

#ifdef USE_SHMAT
//Converts Mat<Type> to Mat<short>
typedef Mat<short> ShMat;
template <class Type> Mat<short> asShMat(const Mat<Type>&);
#endif

#ifdef USE_USHMAT
//Converts Mat<Type> to Mat<unsigned short>
typedef Mat<unsigned short> UShMat;
template <class Type> Mat<unsigned short> asUShMat(const Mat<Type>&);
#endif

#ifdef USE_CHRMAT
//Converts Mat<Type> to Mat<char>
typedef Mat<char> ChrMat;
template <class Type> Mat<char> asChrMat(const Mat<Type>&);
#endif

#ifdef USE_UCHRMAT
//Converts Mat<Type> to Mat<unsigned char>
typedef Mat<unsigned char> UChrMat;
template <class Type> Mat<unsigned char> asUChrMat(const Mat<Type>&);
#endif

#ifdef USE_COMPMAT
typedef Mat<dcomplex> CompMat;
template <class Type> Mat<dcomplex> asCompMat(const Mat<Type>& A);
template <class Type> Mat<dcomplex> asCompMat(const Mat<Type>& Re, const Mat<Type>& Im);
#endif

#ifdef USE_FCOMPMAT
typedef Mat<fcomplex> fCompMat;
template <class Type> Mat<fcomplex> asFcompMat(const Mat<Type>&);
template <class Type> Mat<fcomplex> asFcompMat(const Mat<Type>& Re, const Mat<Type>& Im);
#endif

// Some mathematical operations. These had to be implemented as template
// functions in order to allow different types of matrices to be operated on
// ***************************** Addition ****************************
//
// Some mathematical operations
//
template <class T1, class T2>
Mat<T1>& operator += (Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    std::cerr << "Matrices of incompatible sizes for +=" << std::endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ += (T1) *bPtr++;
  
  return A;
}

template <class T1, class T2>
Mat<T1> operator + (const Mat<T1>& A, const Mat<T2>& B) 
{
  Mat<T1> T(A); 
  return T += B; 
}

// ***************************** Substraction *************************
template <class T1, class T2>
Mat<T1>& operator -= (Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    std::cerr << "Matrices of incompatible sizes for -=" << std::endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ -= (T1) *bPtr++;
  
  return A;
}

template <class T1, class T2>
Mat<T1> operator - (const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return T -= B; }

// ***************************** Multiplication ***********************
template <class T1, class T2>
Mat<T1> operator * (const Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1>& operator *= (Mat<T1>& A, const Mat<T2>& B) {
  A = A * B; return A; }

// ***************************** Division ***** ***********************
template <class T1, class T2>
Mat<T1> operator / (const Mat<T1>& A, const Mat<T2>& B) {
  return A * inv(B); }

template <class T1, class T2>
Mat<T1>& operator /= (Mat<T1>& A, const Mat<T2>& B) {
  A = A * inv(B); return A; }

// *********************** Point multiplication ***********************
template <class T1, class T2>
Mat<T1>& pmultEquals(Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> pmult(const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return pmultEquals(T, B); }

// ************************* Point division ***************************
template <class T1, class T2>
Mat<T1>& pdivEquals(Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> pdiv(const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return pdivEquals(T, B); }

/******************************************************************************
 *     MATRIX CLASS
 ******************************************************************************/

template <class Type> 
class Mat {
  friend class MatrixIterator<Type>;
  friend class ConstMatrixIterator<Type>;

protected:
/*****************************Internal data elements***************************/
   unsigned   _rows;
   unsigned   _cols;
   unsigned   _maxrows;
   unsigned   _maxcols;
   Type	    **_el;
   static unsigned _rangeErrorCount;
/******************************************************************************/
   
public:
  static Boolean flushToDisk;
  enum vector_orientation { ROW = 0, COLUMN = 1}; 

/*******************************Constructors***********************************/
   //Default constructor-try not to use this one. An object instantiated using
   //this constructor cannot be reassigned or copied to
   Mat(void) { _rows = _maxrows = 0; _cols = _maxcols = 0; _el=0; }
   
   //Copy constructor
   Mat(const Mat& A);

  // Create row (default) or column vector from SimpleArray
  // Dir spec removed because DCC can't distinghuish this from Mat(unsigned,unsigned)???
  Mat(const SimpleArray<Type>& A, char dir = Mat<Type>::ROW);
 
   //When the below constructor is used for type character, single quotes must
   // be used for the value which the matrice is to be filled with 
   Mat(unsigned nrows, unsigned ncols, Type value); 

   //Initializes all elements to zero
   Mat(unsigned nrows, unsigned ncols);

   Mat(unsigned nrows, unsigned ncols, const Type *data);
   
   //Initialize object from a matlab file
   Mat(const char *filename, int type = RAW, const char *varname = "A");

/********************************Destructor************************************/
   //Calls the clear function to deallocate memory
  virtual ~Mat();

/***********************Explicit Memory deallocation***************************/
   //May be called within the program to deallocate the dynamic memory
   //assigned to an object.  The object itself is not destroyed until
   //it goes out of scope.
   void clear();

/*******************Element wise accessing operators***************************/
   //Matrix is treated as a linear array with rows being concatenated
   //i.e.  [ 1 2 3     is viewed as [1 2 3 4 5 6]
   //        4 5 6 ]  
   //and single argument n is used to return the appropriate element
   //As with all matrix accessing operators, the first element is element 
   //zero
   //Does allow modification of the element      
   Type& operator() (unsigned n);
   
   //Same as above but does not allow modification of the actual element
   Type  operator() (unsigned n) const;

   //Returns the address of the element of the Matrix specified by row,column
   //Allows modification of the actual element
   Type& operator() (unsigned r, unsigned c);
   
   //Returns the element of the Matrix specified by row,column
   //Does not allow modification of the actual element
   Type  operator() (unsigned r, unsigned c) const;

   //Return a Matrix(Mat) object consisiting of the subsection of the original
   //Matrix(Mat) specified by r1(starting row) to r2(ending row) and c1
   //(starting column) to c2(ending column)
   Mat operator()(unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;
   Mat row(unsigned r) const                { return (*this)(r, r, 0, _cols - 1); }
   Mat rows(unsigned r1, unsigned r2) const { return (*this)(r1, r2, 0, _cols - 1); }
   Mat col(unsigned c) const                { return (*this)(0, _rows - 1, c, c); }
   Mat cols(unsigned c1, unsigned c2) const { return (*this)(0, _rows - 1, c1, c2); }   
/*****************Functions to access Internal data pointers*******************/
   //Returns a double pointer to the first element of the data
   //The primary pointer points to the first of a group of
   //secondary pointers(one for each row).  Each of the secondary pointers
   //points to the first data element for that row
   //The primary pointer that is returned cannot be modified in itself
   //However, a copy of the value of the secondary pointer that it points to
   //can be made and then changes may ne made dirrectly to the encapsulated data
   //This violates the protection provided by the class but may be 
   //necessary in certain places to increase execution speed 
   const Type **getEl() const { return (const Type **) _el; } // Be careful with this one

   //This returns a copy of the secondary pointer to the beginning of the
   //row specified by indx
   const Type* operator[] (unsigned indx) const;

/*********************Exclusive Reassignment Operators*************************/
  //This operator reassigns the calling object to a copy of the argument object
  void  operator() (const Mat&);
  
  //This operator functions identically to the one above execept that it 
  //returns a pointer to the calling object so that assignments maybe cascaded
  Mat& operator = (const Mat& A);
  
  // Functionally identical to operator = (), but directly copies the data
  // pointers from A, thus destroying it (!!). 
  Mat& absorb(Mat& A);

/**************************Matrix reshaping functions**************************/
  // Padding functions; return matrices padded with <value> all around. 
  // Generic pad function; pads to nrows x ncols and inserts at row, col.
  Mat& pad(unsigned nrows, unsigned ncols, int row, int col, Type value = 0);
  Mat  pad(unsigned nrows, unsigned ncols, int row, int col, 
	   Type value = 0) const {
    return padConst(nrows, ncols, row, col, value); }
  Mat  padConst(unsigned nrows, unsigned ncols, int row, int col, 
		Type value = 0) const {
    Mat<Type> A(*this); return A.pad(nrows, ncols, row, col, value); }
  
  // Symmetrical pad, i.e., the size of the returned matrix is (rows+2*rowpad, 
  // cols+2*colpad).
  Mat& pad(unsigned rowpad, unsigned colpad, Type value = 0) {
    return pad(_rows + 2*rowpad, _cols + 2*colpad, rowpad, colpad, value); }
  Mat  pad(unsigned rowpad, unsigned colpad, Type value = 0) const {
    return padConst(rowpad, colpad, value); }
  Mat  padConst(unsigned rowpad, unsigned colpad, Type value = 0) const {
    return padConst(_rows + 2*rowpad, _cols + 2*colpad, rowpad, colpad, value); }
  
  // Resize function. 
  // If the matrix is made smaller, the elements not included in the new matrix
  // will remain in existence (i.e., the memory occupied is not freed).
  // If the matrix is made larger, the existing elements will be preserved.
  // The old matrix elements matrix will be entered at (row, col) of the new
  // resized matrix (row and col may be negative, defaults: 0);
  Mat& resize(unsigned nrows, unsigned ncols, int row = 0, int col = 0);
  
  // This function returns the matrix dimensions to their maximum values
  // before any previous resizing had been done.  Data in the newly resized
  // portion may or may not be the same as what was previously there
  Mat& resize();
  
  //Returns a new matrix which is the side by side concatenation of the calling
  //matrix and the operand
  //No changes are made to the operand or calling matrix object
  Mat appendRight(const Mat& A) const;
  
  //Returns a new matrix which is the top and bottom concatenation of the
  //calling matrix and the operand.
  //No changes are made to the operand or calling matrix object
  Mat appendBelow(const Mat& A) const;
  
  // Swap two rows or columns
  Mat& swapRows(unsigned r1, unsigned r2);
  Mat& swapCols(unsigned c1, unsigned c2);

  //Returns the residual matrix after eliminating row and col
  Mat residual(unsigned row, unsigned col) const;
  
  // Insert a matrix into another at (row, col). If necessary, the argument
  // matrix will be clipped to the bounds of the source matrix.
  Mat& insert(const Mat<Type>& A, int row = 0, int col = 0);
  // Same, but reads the matrix from a raw file
  Mat& insert(const char *path, unsigned nrows=0, unsigned ncols=0, int row=0,int col=0);
  
/*******************Special Matrix reinitialization functions******************/
  // Fills the matrix with <value>
  Mat& fill(Type value);
  
  // Fills a sub-rectangle of the matrix with <value>
  Mat& fill(Type value, unsigned r1, unsigned r2, unsigned c1, unsigned c2);

  // Fills the matrix with <value> in a circular area, centered in the matrix.
  // When <diameter> <= 0, the circle will extend to the edge of the matrix.
  Mat& fillCircle(Type value, double diameter = 0) {
    if (diameter <= 0) diameter = MIN(_rows, _cols);
    return fillEllips(value, diameter, diameter); 
  }

  // Fills the matrix with <value> in the specified circular area, 
  // centered at (row, col). When diameter <= 0, the circle will extend to the
  // closest edge of the matrix.
  Mat& fillCircle(Type value, double row, double col, double diameter = 0) {
    if (diameter <= 0) {
      double rowDiameter = 2*MIN(row + 0.5, _rows - row - 0.5);
      double colDiameter = 2*MIN(col + 0.5, _cols - col - 0.5);
      diameter = MIN(rowDiameter, colDiameter);
    }
    return fillEllips(value, row, col, diameter, diameter); 
  }

  // Fills the matrix with <value> in a circular area, centered in the matrix.
  // When one of the diameters is <= 0, the ellips will extend to the edge of the
  // matrix.
  Mat& fillEllips(Type value, double rowDiameter = 0, double colDiameter = 0);

  // Fills the matrix with <value> in the specified elliptical area, 
  // centered at (row, col). If one of the diameters is <= 0, the ellips will
  // extend to the closest edge of the matrix
  Mat& fillEllips(Type value, double row, double col, double rowDiameter, 
		  double colDiameter);

   //Turns the calling object into a identity matrix
   Mat& eye();
   
   //Fills the calling object with random numbers
   Mat& randuniform(double min = 0, double max = 1);
   
   //Stores a random normal matrix to the calling object
   Mat& randnormal(double mean = 0, double std = 1);
      
   //Stores a hamming window to the calling object
   Mat& hamming();

   //Stores a hanning window to the calling object
   Mat& hanning();

   //Stores a blackman window to the calling object
   Mat& blackman();
   
   //This creates a diagonal matrix from the existing row or column vector objec
   //This function uses(internally) the () operator of this class(The one that 
   //takes one argument and the one that takes two arguments)      
   Mat diag() const;
   
   //This creates a toeplitz matrix using two vectors, the calling object is
   //used as the column vector and the argument as the row vector
   //With no argument, the calling object is used for both
   //This function requires use of the () operator and the appendBelow function
   Mat toeplitz(const Mat& r)const;
   Mat toeplitz() const { return  toeplitz(*this); }

/***********Functions to examine a Matrix(Mat) objects properties**************/

  // Returns non-zero if the matrix is empty
  Boolean operator ! () const { return Boolean(!_rows || !_cols || !_el); }
  operator void * () const    { return (void *) (_rows && _cols && _el); }

   //Returns the number of ACTIVE rows in a Mat object
   unsigned  getrows() const { return _rows; }

   //Returns the number of ACTIVE columns in a Mat object
   unsigned  getcols() const { return _cols; }

   //Returns the number of ACTIVE elements in a Mat object
   unsigned  nElements() const { return _rows*_cols; }

   //Returns the number of ACTUAL rows in a Mat object
   unsigned  getmaxrows() const { return _maxrows; }

   //Returns the number of ACTUAL columns in a Mat object
   unsigned  getmaxcols() const { return _maxcols; }

   //Returns 1 if the ACTIVE Matrix(Mat) is a vector else 0
   int    isvector() const       { return (_rows == 1) || (_cols == 1); }

   //Returns 1 if the ACTIVE Matrix(Mat) is a column vector else 0
   int    iscolumnvector() const { return (_cols == 1); }

   //Returns 1 if the ACTIVE Matrix(Mat) is a row vector else 0
   int    isrowvector() const    { return (_rows == 1); }
  
   //Returns the length of the larger ACTIVE dimension(row vs column) of the
   //Matrix
   unsigned    length() const         { return (_rows > _cols) ? _rows : _cols; }
  
   //Display the contents of the entire ACTIVE Matrix(Mat) object
   std::ostream& display(std::ostream& os) const;
   
   //Display a subsection of the Matrix(Mat) defined by r1 (starting row) to
   //r2(ending row) and c1(starting column) to c2(ending column)
   std::ostream& display(std::ostream& os, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;

/**********************Basic Matrix Arithmetic operators***********************/
  // Unary -
  Mat operator - () const { Mat<Type> A(_rows, _cols, Type(0)); return A -= *this; }

  // Elementwise comparison of two matrices
  int operator != (const Mat&) const;
  int operator == (const Mat& A) const { return !(operator != (A)); }
  
  //Adds a single value to each of the elements of the calling Matrix
  //(on the right hand side) and returns a pointer to this Matrix
  //Right-hand operand is modified
  Mat& operator += (dcomplex);
  
  //Same as above, except that no changes are made to either operand,
  //rather a new Matrix is constructed and returned by the function 
  Mat  operator + (dcomplex x) const      { Mat<Type> A(*this); return A += x; }
  
  //Subtracts a single value from each of the elements of the calling Matrix
  //(on the right hand side) and returns a pointer to this Matrix
  //Right-hand operand is modified
  Mat& operator -= (dcomplex x)            { return *this += (-x); }
  
  //Same as above, except that no changes are made to either operand,
  //rather a new Matrix is constructed and returned by the function
  Mat  operator -  (dcomplex x) const      { return *this +  (-x); }
  
  //Multiplys each element of the right hand side Matrix by a double 
  //on the left hand side.  Changes are propogated to the right hand
  //side Matrix
  Mat& operator *= (dcomplex);
  
  //Same as above, except that changes are not made to either operand
  Mat  operator *  (dcomplex x) const        { Mat<Type> A(*this); return A *= x; }
  
  //Divides each element of the Matrix by a complex and effects the changes in
  //the original Matrix
  Mat& operator /= (dcomplex x)              { return (*this) *=  (1.0/x); }
  
  //Same as above, except that changes are not made to either operand
  Mat  operator /  (dcomplex x) const        { return (*this) *  (1.0/x); }
  
  /********************General purpose Matrix functions**************************/
  //Returns a Matrix which is the transpose of the original matrix
  //Does not change the original matrix
  Mat t() const;
  
  //Returns a Matrix which is the reverse or 180 rotate 
  //Does not change the original matrix
  Mat rotate180() const;
  
  //Returns a Matrix which is the inverse of the original matrix
  //Does not change the original matrix. Returns an empty matrix on failure
  Mat inv() const;
  
  //Returns the hermitan transpose of the original Matrix
  //Does not change the original matrix
  Mat h() const;
  
  //Returns the value of the smallest element of the ACTIVE matrix
  //optionally returns the index (row and column) of the element
  Type min(unsigned *row = 0, unsigned *col = 0) const;
  
  //Returns the value of the largest element of the ACTIVE matrix
  //optionally returns the index (row and column) of the element
  Type max(unsigned *row = 0, unsigned *col = 0) const;
  
  // Returns the median of the matrix, only considering elements in the supplied range.
  // By default, the full range of the matrix will be used.
  Type median(Type minVal = 0, Type maxVal = 0) const;

  //Returns the sum of all the elements of the ACTIVE matrix
  double  sum() const { return real(csum()); }
  dcomplex csum() const;
  
  //Returns the sum of squares of all the elements of the ACTIVE matrix
  double  sum2() const { return real(csum2()); }
  dcomplex csum2() const;
  
  //Returns the sum of all the elements of the ACTIVE matrix divided by the
  //number of active elements
  double  mean() const  { return sum()/nElements(); }
  dcomplex cmean() const { return csum() / double(nElements()); }

  //Returns the standard deviation of the matrix elements
  double  std() const  { return ::sqrt(var()); }
  dcomplex cstd() const { return std::sqrt(cvar()); }

  //Returns the variance of the matrix elements
  double  var() const  { double  mn = mean(); return sum2()/nElements() - SQR(mn); }
  dcomplex cvar() const { dcomplex mn = cmean(); return csum2() / double(nElements()) - SQR(mn); }
  
  //Returns the norm of the matrix
  double  norm() const  { return ::sqrt(sum2()); }
  dcomplex cnorm() const { return std::sqrt(csum2()); }
  
  //Returns the trace of the matrix
  double  trace() const { return real(ctrace()); }
  dcomplex ctrace() const;
  
  //Returns the determinant of the matrix
  double  det() const { return real(cdet()); }
  dcomplex cdet() const;
  
  /**********************Other Matrix element wise functions*********************/
  //Applies a given single argument function to each of the elements
  Mat&  applyElementWise(double (*function)(double));
  Mat   applyElementWise(double (*function)(double)) const {
    return applyElementWiseConst(function); }
  Mat   applyElementWiseConst(double (*function)(double)) const {
    Mat<Type> A(*this); return A.applyElementWise(function); }

  Mat&  applyElementWiseC2D(double (*function)(dcomplex));
  Mat   applyElementWiseC2D(double (*function)(dcomplex)) const {
    return applyElementWiseConstC2D(function); }
  Mat   applyElementWiseConstC2D(double (*function)(dcomplex)) const {
    Mat<Type> A(*this); return A.applyElementWiseC2D(function); }

  Mat&  applyElementWiseC2C(dcomplex (*function)(dcomplex));
  Mat   applyElementWiseC2C(dcomplex (*function)(dcomplex)) const {
    return applyElementWiseConstC2C(function); }
  Mat   applyElementWiseConstC2C(dcomplex (*function)(dcomplex)) const {
    Mat<Type> A(*this); return A.applyElementWiseC2C(function); }
  
  Mat& applyIndexFunction(IndexFunction F);
  Mat  applyIndexFunction(IndexFunction F) const {
    return applyIndexFunctionConst(F); }
  Mat   applyIndexFunctionConst(IndexFunction F) const {
    Mat<Type> A(*this); return A.applyIndexFunction(F); }

  Mat& applyIndexFunction(ComplexIndexFunction F);
  Mat  applyIndexFunction(ComplexIndexFunction F) const {
    return applyIndexFunctionConst(F); }
  Mat   applyIndexFunctionConst(ComplexIndexFunction F) const {
    Mat<Type> A(*this); return A.applyIndexFunction(F); }

  //Takes the inverse natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& exp();
  Mat  exp() const      { return expConst(); }
  Mat  expConst() const { Mat<Type> T(*this); return T.exp(); }
  
  //Takes the natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& log();
  Mat  log() const      { return logConst(); }
  Mat  logConst() const { Mat<Type> T(*this); return T.log(); }
  
  //Takes the cosine of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& cos();
  Mat  cos() const      { return cosConst(); }
  Mat  cosConst() const { Mat<Type> T(*this); return T.cos(); }
  
  //Takes the sine of each of the elements of the matrix
  Mat& sin();
  Mat  sin() const      { return sinConst(); }
  Mat  sinConst() const { Mat<Type> T(*this); return T.sin(); }
  
  //Returns the conjugate matrix
  Mat& conj();
  Mat  conj() const      { return conjConst(); }
  Mat  conjConst() const { Mat<Type> T(*this); return T.conj(); }
  
  //Takes the absolute value of each of the elements of the matrix
  Mat& abs();
  Mat  abs() const      { return absConst(); }
  Mat  absConst() const { Mat<Type> T(*this); return T.abs(); }
  
  //Rounds each of the elements of the matrix
  Mat& round();
  Mat  round() const      { return roundConst(); }
  Mat  roundConst() const { Mat<Type> T(*this); return T.round(); }

  //Takes the sqrt value of each of the elements of the matrix
  Mat& sqrt();
  Mat  sqrt() const      { return sqrtConst(); }
  Mat  sqrtConst() const { Mat<Type> T(*this); return T.sqrt(); }

  //Takes the power value of each of the elements of the matrix
  Mat& pow(double exponent);
  Mat  pow(double exponent) const      { return powConst(exponent); }
  Mat  powConst(double exponent) const { Mat<Type> T(*this); return T.pow(exponent); }

/************************Special purpose Matrix functions**********************/
  // Returns transpose(A)*A of the matrix. 
  Mat transposeXself() const;

   void   eig(Mat& D, Mat& V) const;

   //Returns a Matrix which is the householder transform of the calling object
   //Calling object must be a vector
   Mat house() const;

   //Returns a Matrix which is the rowhouse transform of the calling object and
   //the input column vector
   Mat rowhouse(const Mat& V) const;

   //convolve 2d Matrix //taken from red book and modified
   Mat convolv2d(const Mat& filter) const;

  //Returns a histogram for the calling object using the specified range and # bins
  Histogram histogram(double minin = 0, double maxin = 0, unsigned n = 0) const;

  Mat& histmod(const Histogram& hist1, const Histogram& hist2);
  Mat  histmod(const Histogram& hist1, const Histogram& hist2) const {
    return histmodConst(hist1, hist2); }
  Mat  histmodConst(const Histogram& hist1, const Histogram& hist2) const {
    Mat<Type> A(*this); return A.histmod(hist1, hist2); }

  // Returns an array which contains only the values in the range specified
  SimpleArray<Type> array(Type minVal = 0, Type maxVal = 0) const;

  // Returns first (topleft) element
  Type              scalar() const { return **_el; }

   //Performs a QR factorization and modifies Q and R appropriately
   //In this case, X=Q*R
   //Where X is the calling object
   void   qr(Mat& R, Mat& Q) const;

   //Performs a QR factorization using pivoting.
   //In this case X*P=Q*R
   //Where X is the calling object
   void   qr(Mat& R, Mat& Q, Mat& P) const;

   //Performs the singular value decomposition on the calling matrix
   //and modifies U, S, and V appropriately
   //In this case X = U*S*(V')
   //Where X is the calling object
   void   svd(Mat& U, Mat& S, Mat& V) const;

  // Sets all values in the matrix below minVal to minFill, and those
  // above maxVal to maxFill.
  Mat& clip(Type minVal, Type maxVal, Type minFill, Type maxFill);
  Mat  clip(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    return clipConst(minVal, maxVal, minFill, maxFill); }
  Mat  clipConst(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    Mat<Type> A(*this); return A.clip(minVal, maxVal, minFill, maxFill); }

  // Maps the values of a matrix through a value map
  Mat& map(const ValueMap& valueMap);
  Mat  map(const ValueMap& valueMap) const { return mapConst(valueMap); }
  Mat  mapConst(const ValueMap& valueMap) const { 
    Mat<Type> A(*this); return A.map(valueMap);}

  // Scales a Matrix (similar to map, but truncates the output range to be
  // [minout, maxout])
  Mat& scale(double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin = 0.0);
  Mat  scale(double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin = 0.0) const {
    return scaleConst(minout, maxout, minin, maxin); }
  Mat  scaleConst(double minout = 0.0, double maxout = 255.0, 
		  double minin = 0.0, double maxin = 0.0) const {
    Mat<Type> A(*this); return A.scale(minout, maxout, minin, maxin); }

   //Performs a linear interpolation on the calling object Matrix
   //The y (argument) Matrix recieves the results
   //Urow is the upscale factor for the rows and Drow is the downscale factor
   //for the rows.  Ucol and Dcol are analagous for the columns.
   void linearinterpolate( unsigned Urow, unsigned Drow, unsigned Ucol, unsigned Dcol,Mat& y) const;

   //Applies the filter defined by B and A to the calling object and returns the
   //resultant Matrix.
   //The filter has the form y= (B/A)*x
   //Where x is the calling object and y is the result
   Mat filter(Mat& B,Mat& A) const;

/***************************morphology functions*******************************/
   // In the following (greylevel) morphology operations, negative strel
   // values are considered NaN's and ignored.
#ifdef USE_DBLMAT
   Mat erode(const Mat<double>& strel) const;
   Mat dilate(const Mat<double>& strel) const;
   Mat open(const Mat<double>& strel) const  { return erode(strel).dilate(strel); }
   Mat close(const Mat<double>& strel) const { return dilate(strel).erode(strel); }
#endif

  //   Not converted yet
  //   void radialscale(Mat& scale, Mat& y) const;

  //This function converts a Matlab matrix object to a Mat object
  //The first argument is a pointer to a Matlab matrix object
  //This function checks the storage format of the Matlab object
  //(Full or Sparse) and does the appropriate conversions to return a Mat
  //object of the calling object's type.
#ifdef HAVE_MATLAB
  Mat matlab2Mat(Matrix *MatlabMatrix) const;
  
  //Creates a matlab type matrix and assigns it the name specified as the
  //argument.  It then returns a pointer to this object
  Matrix *mat2Matlab(const char *name) const;
#endif

/*************************Matrix file I/O functions****************************/
   //This function goes to the named mat file and searches for the given
   //variable.  Next, if possible, it reinitializes the calling object to
   //represent the same matrix as that defined by the Matlab variable.  The
   //Matlab matrix is stored in double format.  The function Matlab2Mat will
   //convert the Matlab data to the type of its calling object (this).  This may
   //mean a loss of precision.
#ifdef HAVE_MATLAB
  Boolean loadMatlab(const char *filename, const char *varname = "A");
#endif

   //Loads a binary version of a matrix from a file
   //must have previously been stroed by saveRaw
   //The calling object must have the same type as the object that was 
   //previously saved in order for this function to work correctly
  Boolean loadRaw(const char *filename, unsigned nrows = 0, unsigned ncols = 0);

   //Ascii files will always have rows and columns up top
   //Loads an Ascii version of the saved Matrix
  Boolean loadAscii(const char *filename);
   
   //Loads a matrix from either a binary, Ascii or Matlab file, determined
   //by the type flag
  Boolean load(const char *filename, int type = RAW, const char *varname = "A");
   //Saves the calling object in either binary,Ascii or Matlab style, determined
   //by a flag set in this header file
  Boolean save(const char *filename, int type = RAW, const char *varname = "A") const;

   //Saves the calling object as a Matlab readable matrix to the namedfile with
   //given variable name.  Matlab stores in double format so, the calling object
   //is first converted to a temporary Mat object of type double.
   //Matlab stores in a column based format, so the Mat<double> is transposed
   //before saving.  The third argument lets the user update a file("u") or
   //rewrite an existing file or create a new one("w") 
#ifdef HAVE_MATLAB
  Boolean saveMatlab(const char *filename, const char *varname = "A", 
		     const char *option = "u") const;
#endif

   //Saves the calling matrix in binary format
  Boolean saveRaw(const char *) const;

   //Saves the calling object in Ascii format
   //Ascii files will always save rows and columns up top 
  Boolean saveAscii(const char *) const;
  
/*****************************fft functions************************************/

  // For the following fft functions, a dimension of 0 (which is the
  // default argument) will be expanded to the next power of two, equal
  // to or above the matrix dimension. A dimension of 1 indicates that
  // the fft should not be taken along that dimension.

  // Non-const fft function
  Mat& fft(unsigned nrows = 0, unsigned ncols = 0) {
    return _fft(nrows, ncols, ::fft); }
  // Const fft function
  Mat  fft(unsigned nrows = 0, unsigned ncols = 0) const { 
    return fftConst(nrows, ncols); }
  // Explicit const fft function
  Mat  fftConst(unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat<Type> A(*this); return A.fft(nrows, ncols); }

  // Non-const ifft function
  Mat& ifft(unsigned nrows = 0, unsigned ncols = 0);
  // Const ifft function
  Mat  ifft(unsigned nrows = 0, unsigned ncols = 0) const { 
    return ifftConst(nrows, ncols); }
  // Explicit const ifft function
  Mat  ifftConst(unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat<Type> A(*this); return A.ifft(nrows, ncols); }

/****************************End of public functions***************************/

private:

/********************************Private functions*****************************/
   //Allocate the memory for a Matrix object
   //If the object already has memory allocated to it and this function is 
   //called, it first deallocates the memory of the object
   void _allocateEl();

   //This function copies the Matrix argument B to the subsection of
   //the calling object defined by r1(starting row),r2(ending row),
   //c1(starting column),c2(ending column).  
   void _modify_sub_section(unsigned r1,unsigned r2,unsigned c1,unsigned c2,const Mat& B);

  // Given the path of a raw file, attempts to infer the matrix dimensions from
  // the size of the file.
  void _checkMatrixDimensions(const char *path, unsigned& nrows, unsigned& ncols) const;

  Mat& _fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc);

/***************************End Private functions*****************************/
   
/*******************************END OF Mat CLASS******************************/
};

// Non-member template functions

template <class T>
void clear(Mat<T> A) { A.clear(); }

template <class T>  
Mat<T> pad(const Mat<T>& A, unsigned rowpad, unsigned colpad, T value = 0)
{
  return A.padConst(rowpad, colpad, value);
}

template <class T>  
Mat<T>& resize(Mat<T> A, unsigned nrows, unsigned ncols, 
	       int row = 0, int col = 0)
{
  return A.resize(nrows, ncols, row, col);
}

template <class T>
Mat<T>& resize(Mat<T> A) { return A.resize(); }

template <class T>  
Mat<T> appendRight(const Mat<T>& A, const Mat<T>& V)
{
  return A.appendRight(V);
}

template <class T>
Mat<T> appendBelow(const Mat<T>& A, const Mat<T>& V)
{
  return A.appendBelow(V);
}
  
template <class T> 
Mat<T> residual(const Mat<T>& A, unsigned row, unsigned col)
{
  return A.residual(row, col);
}
  
template <class T>
Mat<T>& fill(Mat<T> A, T value) { return A.fill(value); }

template <class T>
Mat<T>& eye(Mat<T>& A) { return A.eye(); }

template <class T>
Mat<T>& randuniform(Mat<T>& A, double min = 0, double max = 1)
{
  return A.randuniform(min, max);
}

template <class T>
Mat<T>& randnormal(Mat<T>& A, double mean = 0, double std = 1)
{ 
  return A.randnormal(mean, std);
}

template <class T>
Mat<T>& hamming(Mat<T>& A) { return A.hamming(); }

template <class T>
Mat<T>& hanning(Mat<T>& A) { return A.hanning(); }

template <class T>
Mat<T>& blackman(Mat<T>& A) { return A.blackman(); }

template <class T>
Mat<T> diag(const Mat<T>& A) { return A.diag(); }

template <class T>
Mat<T> toeplitz(const Mat<T>& c, const Mat<T>& r)
{ 
  return(c.toeplitz(r));
}

template <class T>
Mat<T> toeplitz(const Mat<T>& c) { return(c.toeplitz(c)); }

template <class T>
unsigned getrows(const Mat<T>& A) { return A.getrows(); }

template <class T>
unsigned getcols(const Mat<T>& A) { return A.getcols(); }

template <class T>
unsigned nElements(const Mat<T>& A) { return A.nElements(); }

template <class T>
unsigned getmaxrows(const Mat<T>& A) { return A.getmaxrows(); }

template <class T>
unsigned get_maxcols(const Mat<T>& A) { return A.getmaxcols(); }

template <class T>
int isvector(const Mat<T>& A) { return A.isvector(); }

template <class T>
int iscolumnvector(const Mat<T>& A) { return A.iscolumnvector(); }

template <class T>
int isrowvector(const Mat<T>& A) { return A.isrowvector(); }

template <class T>
unsigned length(const Mat<T>& A) { return A.length(); }

template <class T>
std::ostream& display(std::ostream& os, const Mat<T>& A) { return A.display(os); }

template <class T>
std::ostream& display(std::ostream& os, const Mat<T> A,unsigned r1, unsigned r2, 
		 unsigned c1, unsigned c2)
{ 
  return A.display(os, r1, r2, c1, c2); 
}
  
template <class T>
Mat<T>  t(const Mat<T>& A) { return A.t(); }

template <class T>
Mat<T>  rotate180(const Mat<T>& A) { return A.rotate180(); }

template <class T>
Mat<T>  h(const Mat<T>& A) { return A.h(); }
  
template <class T>
double sum(const Mat<T>& A) { return A.sum(); }

template <class T>
dcomplex csum(const Mat<T>& A) { return A.csum(); }

template <class T>
double sum2(const Mat<T>& A) { return A.sum2(); }

template <class T>
dcomplex csum2(const Mat<T>& A) { return A.csum2(); }

template <class T>
double mean(const Mat<T>& A) { return A.mean(); }

template <class T>
dcomplex cmean(const Mat<T>& A) { return A.cmean(); }

template <class T>
double trace(const Mat<T>& A) { return A.trace(); }

template <class T>
dcomplex ctrace(const Mat<T>& A) { return A.ctrace(); }

template <class T>
double norm(const Mat<T>& A) { return A.norm(); }

template <class T>
dcomplex cnorm(const Mat<T>& A) { return A.cnorm(); }

template <class T>
double det(const Mat<T>& A) { return A.det(); }

template <class T>
dcomplex cdet(const Mat<T>& A) { return A.cdet(); }

template <class T>
Mat<T> exp(const Mat<T>& A) { return A.exp(); }

template <class T>
Mat<T> sqrt(const Mat<T>& A) { return A.sqrt(); }

template <class T>
Mat<T> rowhouse(const Mat<T>& A, const Mat<T>& V)
{
  return A.rowhouse(V);
}

template <class T>
Mat<T> convolv2d(const Mat<T>& A, const Mat<T>& filter) 
{
  return A.convolv2d(filter);
}
  
template <class T>
void remake(Mat<T>& A, unsigned nrows, unsigned ncols)
{
  A(Mat<T>(nrows, ncols));
}
  
#ifdef USE_COMPMAT
#ifdef USE_DBLMAT
extern Mat<double> applyElementWiseC2D(const Mat<dcomplex>& A,
				       double (*function)(const dcomplex&));
extern Mat<double> arg(const Mat<dcomplex>& A);
extern Mat<double> real(const Mat<dcomplex>& A);
extern Mat<double> imag(const Mat<dcomplex>& A);
#endif
#endif
#ifdef USE_FCOMPMAT
#ifdef USE_FLMAT
extern Mat<float> applyElementWiseC2D(const Mat<fcomplex>& A,
					double (*function)(const dcomplex&));
extern Mat<float> arg(const Mat<fcomplex>& A);
extern Mat<float> real(const Mat<fcomplex>& A);
extern Mat<float> imag(const Mat<fcomplex>& A);
#endif
#endif

template <class T> 
Mat<T>& operator += (double addend, Mat<T>& A) { return A+=addend; }

template <class T> 
Mat<T>  operator +  (double addend, const Mat<T>& A) { return A+addend; }

template <class T> 
Mat<T>& operator -= (double subend, Mat<T>& A) { return A-=subend; }

template <class T> 
Mat<T>  operator -  (double subend, const Mat<T>& A) { return A-subend; }

template <class T> 
Mat<T>& operator *= (double factor,  Mat<T>& A) { return A*=factor; }

template <class T> 
Mat<T>  operator *  (double factor,  const Mat<T>& A) { return A*factor; }

template <class T> 
Mat<T>& operator /= (double factor,  Mat<T>& A) { return A/=factor; }

template <class T> 
Mat<T>  operator /  (double factor,  const Mat<T>& A) { return A/factor; }

template <class T> 
Mat<T> inv(const Mat<T>& A) { return A.inv(); }

template <class T> 
T min(const Mat<T>& A, unsigned *row = 0, unsigned *col = 0)
{ 
  return A.min(row, col);
}

template <class T> 
T max(const Mat<T>& A, unsigned *row = 0, unsigned *col = 0) 
{ 
  return A.max(row, col);
}

template <class T> T 
median(const Mat<T>& A, T minVal, T maxVal) 
{ 
  return A.median(minVal, maxVal);
}

template <class T> 
Mat<T> applyElementWise(const Mat<T>& A, double (*function)(double))
{ 
  return A.applyElementWise(function);
}

template <class T> 
Mat<T> applyElementWiseC2C(const Mat<T>& A, dcomplex (*function)(dcomplex))
{
  return A.applyElementWiseC2C(function);
}

template <class T> 
Mat<T> log(const Mat<T>& A) { return A.log(); }

template <class T> 
Mat<T> cos(const Mat<T>& A) { return A.cos(); }

template <class T>
Mat<T> sin(const Mat<T>& A) { return A.sin(); }

template <class T>
Mat<T> round(const Mat<T>& A) { return A.round(); }

template <class T>
Mat<T> abs(const Mat<T>& A) { return A.abs(); }

template <class T> 
Mat<T> pow(const Mat<T>& A, double exp) { return A.pow(exp); }

template <class T>
void eig(const Mat<T>& A, Mat<T>& D, Mat<T>& V) { A.eig(D, V); }

template <class T> 
Mat<T> house(const Mat<T>& A) { return A.house(); }

template <class T> 
Histogram histogram(const Mat<T>& A, double minin = 0, double maxin = 0, unsigned n = 0) { return A.histogram(minin, maxin, n); }

template <class T>
void qr(const Mat<T>& A, Mat<T>& R, Mat<T>& Q) { A.qr(R, Q); }

template <class T>
void svd(const Mat<T>& A, Mat<T>& U, Mat<T>& S, Mat<T>& V) { A.svd(U,S,V); }
  
template <class T>
Mat<T> clip(const Mat<T>& A, T minVal, T maxVal, T minFill, T maxFill) 
{
  return A.clipConst(minVal, maxVal, minFill, maxFill);
}

template <class T>
Mat<T> map(const Mat<T>& A, const ValueMap& valueMap) 
{ 
  return A.mapConst(valueMap); 
}

template <class T>
Mat<T> scale(const Mat<T>& A, double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin=0.0) 
{ 
  return A.scaleConst(minout, maxout, minin, maxin);
}
  
template <class T>
void linearinterpolate(const Mat<T>& x, unsigned Urow, unsigned Drow, 
		       unsigned Ucol, unsigned Dcol, Mat<T>& y)
{
  x.linearinterpolate(Urow,Drow,Ucol,Dcol,y); 
}
  
template <class T>
Mat<T> filter(Mat<T>& B, Mat<T>& A, const Mat<T>& x) 
{ 
  return (x.filter(B, A));
}
  
#ifdef USE_DBLMAT
template <class T>
Mat<T> erode(const Mat<T>& A, const Mat<double>& strel) 
{
  return A.erode(strel);
}

template <class T>
Mat<T> dilate(const Mat<T>& A, const Mat<double>& strel) 
{
  return A.dilate(strel);
}

template <class T> 
Mat<T> open(const Mat<T>& A, const Mat<double>& strel)
{
  return A.open(strel);
}

template <class T> 
Mat<T> close(const Mat<T>& A, const Mat<double>& strel)
{
  return A.close(strel);
}
#endif // USE_DBLMAT

template <class T>
T scalar(const Mat<T>& A) { return A.scalar(); }
  
template <class T> 
void qr(const Mat<T>& A, Mat<T>& R, Mat<T>& Q, Mat<T>& P) { A.qr(R,Q,P); }
  
  
template <class T>
Boolean loadRaw(Mat<T>& A, const char *filename,
		unsigned nrows = 0, unsigned ncols = 0) 
{
  return A.loadRaw(filename, nrows, ncols);
}

template <class T>
Boolean loadAscii(Mat<T>& A, const char *filename)
{ 
  return A.loadAscii(filename);
}
  
#ifdef HAVE_MATLAB   
template <class T>
Boolean loadMatlab(Mat<T>& A, const char *filename, const char *varname = "A")
{
  return A.loadMatlab(filename, varname);
}
  
template <class T>
Boolean saveMatlab(const Mat<T>& A, const char *filename, 
		   const char *varname, const char *option = "u")
{
  return A.saveMatlab(filename, varname, option);
}
#endif // HAVE_MATLAB
  
template <class T>
Boolean saveRaw(const Mat<T>& A, const char *filename) 
{ 
  return A.saveRaw(filename);
}

template <class T>
Boolean saveAscii(const Mat<T>& A, const char *filename)
{
  return A.saveAscii(filename);
}

#ifdef USE_COMPMAT
template <class T>
Mat<dcomplex> fft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
template <class T>
Mat<dcomplex> ifft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template <class T>
Mat<fcomplex> ffft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
template <class T>
Mat<fcomplex> fifft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
#endif /* USE_FCOMPMAT */

template <class T> std::ostream& operator << (std::ostream&, const Mat<T>&);

/*******************************Sub classes************************************/

template <class Type>
class Ones : public Mat<Type> {
public: 
  Ones(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols, (Type) 1) {}
  Ones(const Mat<Type>& A) : Mat<Type>(A.getrows(), A.getcols(), (Type) 1) {}
  Ones(unsigned n) : Mat<Type>(n, n, (Type) 1) {}
  ~Ones() {}
};


template <class Type>
class Zeros : public Mat<Type> {
public:
  Zeros(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols) {}
  Zeros(const Mat<Type>& A) : Mat<Type>(A.getrows(), A.getcols()) {}  
  Zeros(unsigned n)  : Mat<Type>(n, n) {} 
  ~Zeros() {}
}; 


template <class Type>
class Eye : public Mat<Type> {
public:
  Eye(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols) { this->eye(); }
  Eye(Mat<Type>& A) :  Mat<Type>(A.getrows(), A.getcols()) { this->eye(); }
  Eye(unsigned n) : Mat<Type>(n, n)  { this->eye(); }
  ~Eye() {}
};


template <class Type>
class Randuniform : public Mat<Type> {
public:
  Randuniform(unsigned nrows, unsigned ncols, double min = 0, double max = 1) 
    : Mat<Type>(nrows, ncols) { this->randuniform(min, max); }
  Randuniform(const Mat<Type>& A, double min = 0, double max = 1) 
    : Mat<Type>(A.getrows(), A.getcols()) { this->randuniform(min, max); }  
  Randuniform(unsigned n, double min = 0, double max = 1) : Mat<Type>(n, n) { 
      this->randuniform(min, max); } 
  ~Randuniform() {}
};


template <class Type>
class Randnormal : public Mat<Type> {
public:
  Randnormal(unsigned nrows, unsigned ncols, double mean = 0, double std = 1) 
    : Mat<Type>(nrows, ncols) { this->randnormal(mean, std); }
  Randnormal(const Mat<Type>& A, double mean = 0, double std = 1) 
    : Mat<Type>(A.getrows(), A.getcols()) { this->randnormal(mean, std); }  
  Randnormal(unsigned n, double mean = 0, double std = 1)  
    : Mat<Type>(n, n) { this->randnormal(mean, std); } 
  ~Randnormal() {}
};

template <class Type>
class  Hamming: public Mat<Type> {
public:
  Hamming(unsigned n)  : Mat<Type>(n, 1) { this->hamming(); } 
  ~Hamming() {}
};

template <class Type>
class  Hanning: public Mat<Type> {
public:
   Hanning(unsigned n)  : Mat<Type>(n, 1) { this->hanning(); } 
  ~Hanning() {}
};

template <class Type>
class  Blackman: public Mat<Type> {
public:
   Blackman(unsigned n)  : Mat<Type>(n, 1) { this->blackman(); } 
  ~Blackman() {}
};


/******************************************************************************
 * MATRIX ITERATOR CLASSES
 *
 * Quick traversal of a matrix. 
 * IMPORTANT NOTE: For efficiency reasons, no range cheking is performed; I.e.,
 * the user is responsible for staying within the matrix.
 * 
 ******************************************************************************/

template<class Type>
class MatrixIterator {
  Mat<Type>& _mat;   // matrix to be iterated
  Type      *_elPtr; // Current location in matrix

public:
  MatrixIterator(Mat<Type>& mat) : _mat(mat) { first(); }

  // Reset to first element
  Type& first() { return reset(); }
  // Reset to last element
  Type& last()  { return reset(this->_rows - 1, this->_cols - 1); } 
  // Reset to element (r, c)
  Type& reset(unsigned r=0, unsigned c=0) { return *(_elPtr = (_mat._el)[r] + c); } 

  Type& operator ++()    { return *++_elPtr; } // Prefix increment
  Type& operator ++(int) { return *_elPtr++; } // Postfix increment
  Type& operator --()    { return *--_elPtr; } // Prefix decrement
  Type& operator --(int) { return *_elPtr--; } // Postfix decrement
};

template<class Type>
class ConstMatrixIterator {
  const Mat<Type>& _mat;   // matrix to be iterated
  const Type      *_elPtr; // Current location in matrix

public:
  ConstMatrixIterator(const Mat<Type>& mat) : _mat(mat) { first(); }

  // Reset to first element
  const Type& first() { return reset(); }
  // Reset to last element
  const Type& last()  { return reset(this->_rows - 1, this->_cols - 1); } 
  // Reset to element (r, c)
  const Type& reset(unsigned r=0, unsigned c=0) { return *(_elPtr = (_mat._el)[r] + c); }

  const Type& operator ++()    { return *++_elPtr; } // Prefix increment
  const Type& operator ++(int) { return *_elPtr++; } // Postfix increment
  const Type& operator --()    { return *--_elPtr; } // Prefix decrement
  const Type& operator --(int) { return *_elPtr--; } // Postfix decrement
};

template<class Type>
SimpleArray<Type> array(const Mat<Type>& A, Type minVal = 0, Type maxVal = 0);


//VF: fixing stupid habit of putting templates into .cc file
//
// Some mathematical operations
//

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

template <class T1, class T2>
Mat<T1>& pmultEquals(Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    cerr << "Matrices of incompatible sizes for pmultEquals" << endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ *= (T1) *bPtr++;
  
  return A;
}

template <class T1, class T2>
Mat<T1>& pdivEquals(Mat<T1>& A, const Mat<T2>& B)
{
  unsigned nrows = A.getrows();
  unsigned ncols = A.getcols();

  if (!(A.isvector() && B.isvector() && (A.length() == B.length())) &&
      ((B.getrows() != nrows) || (B.getcols() != ncols))) {
    cerr << "Matrices of incompatible sizes for pdivEquals" << endl;
    return A;
  }

  T1       *aPtr = (T1 *) A.getEl()[0];
  const T2 *bPtr = B.getEl()[0];

  for (unsigned i = nrows; i; i--)
    for (unsigned j = ncols; j; j--)
      *aPtr++ /= (T1) *bPtr++;
  
  return A;
}

/*********************************
Mat class definitions and member functions
**********************************/
//
//-------------------------// 
//

template <class Type>
Mat<Type>::Mat(const Mat<Type>& A)
{
  _rows = A._rows; _cols = A._cols;
  _maxrows = A._maxrows; _maxcols = A._maxcols;
  _el = 0;
   
  _allocateEl();

  if (_maxrows && _maxcols && _el)
    memcpy(*_el, *A._el, _maxrows*_maxcols*sizeof(Type));
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>::Mat(const SimpleArray<Type>& A, char dir)
{
  unsigned N = A.size();
  
  if (dir == Mat<Type>::ROW) {
    _rows = _maxrows = 1;
    _cols = _maxcols = N;
    _el = 0;
    _allocateEl();
  }
  else {
    _rows = _maxrows = N;
    _cols = _maxcols = 1;
    _el = 0;
    _allocateEl();
  }
  
  Type *elPtr = *_el;
  for (unsigned i = 0; i < N; i++)
    *elPtr++ = A[i];
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols, Type value)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
   
   if (value != 0.0)
      fill(value);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(unsigned nrows, unsigned ncols)
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
Mat<Type>::Mat(unsigned nrows, unsigned ncols, const Type *data)
{
   _rows = _maxrows = nrows;
   _cols = _maxcols = ncols;
   _el = 0;
   
   _allocateEl();
   
   size_t  nBytesInRow = ncols*sizeof(Type);
   const Type *dataPtr = data;
   for (unsigned row = 0; row < nrows; row++) {
      memcpy(_el[row], dataPtr, nBytesInRow);
      dataPtr += ncols;
   }
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>::Mat(const char *filename, int type, const char *varname)
{
  _rows = _maxrows = 0;
  _cols = _maxcols = 0;
  _el = 0;
  
  load(filename, type, varname);
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
      if( _el[0] ) {
        delete [] _el[0];
      }
      delete [] _el; // delete row of pointers  
      _el = 0;
   }

   _rows = _maxrows = _cols = _maxcols = 0;
}

/****************************************************************************
These overloaded operators allow accessing of individual elements within the
matrix.  By convention, the first element of the matrix is labeled element 
zero.  Since the matrix is stored in row major format, The element numbers increase going to the right and at the following row. The element at the
last row and last column corresponds to the largest possible value of the 
input argument.
*****************************************************************************/


/****************************************************************************
This accessing operator allows modification of individual elements of the 
matrix. A seperate function is used for accessing elements within
complex matrices.
*****************************************************************************/
template <class Type>
Type& 
Mat<Type>::operator () (unsigned n) 
{
  if (n >= _rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _rows*_cols - 1;
  }

  unsigned int rowindex;
  unsigned int colindex;
  rowindex= (n/_cols);
  colindex= (n%_cols);
  return _el[rowindex][colindex];
}

/****************************************************************************
This accessing operator prevents accidental modification by having a constant
output argument. A seperate function is used for accessing elements within
complex matrices.
*****************************************************************************/

template <class Type>
Type
Mat<Type>::operator () (unsigned n) const 
{
  if (n >= _rows*_cols) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << n << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows*_cols - 1 << endl;
      _rangeErrorCount--;
    }
    n = _rows*_cols - 1;
  }

  unsigned int rowindex;
  unsigned int colindex;
  rowindex= (n/_cols);
  colindex= (n%_cols);
  return _el[rowindex][colindex];
}

//
//-------------------------// 
//
template <class Type>
Type& 
Mat<Type>::operator () (unsigned r, unsigned c) {
  using std::min;		// (bert) force it to find the right min()

  if ((r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << r << ", " << c << ") exceed matrix dimensions. "
	   << "Changed to (" << min(r, _rows - 1) << ", " << min(c, _cols - 1)
	   << ")" << endl;
      _rangeErrorCount--;
    }
    r = min(r, _rows - 1);
    c = min(c, _cols - 1);
  }

  return(_el[r][c]);
}

//
// 
//
template <class Type>
Type
Mat<Type>::operator () (unsigned r, unsigned c) const {
  using std::min;		// (bert) force it to find the right min().

  if ((r >= _rows) || (c >= _cols)) {
    if (_rangeErrorCount) {
      cerr << "Error: indices (" << r << ", " << c << ") exceed matrix dimensions. "
	   << "Changed to (" << min(r, _rows - 1) << ", " << min(c, _cols - 1)
	   << ")" << endl;
      _rangeErrorCount--;
    }
    r = min(r, _rows - 1);
    c = min(c, _cols - 1);
  }

  return(_el[r][c]);
}

//
//------------------// 
//
//cropt a section of the matrix
template <class Type>
Mat<Type>
Mat<Type>::operator () (unsigned r1, unsigned r2, unsigned c1, unsigned c2) const {
  if ((r1 > r2) || (c1 > c2) || (r2 >= _rows) || (c2>= _cols)){
    cerr << "Error in cropping: improper row or column sizes." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }
  Mat<Type> A(r2-r1+1,c2-c1+1);
  Type *Aptr = A._el[0];
  for (unsigned i=r1; i<=r2 ; i++){
    for (unsigned j=c1; j<=c2; j++)
      *Aptr++=_el[i][j];
  }
  return(A);     
}

//
//-------------------------// 
//
template <class Type>
const Type * 
Mat<Type>::operator [] (unsigned index) const {
  if (index >= _rows) {
    if (_rangeErrorCount) {
      cerr << "Error: index " << index << " exceeds matrix dimensions. ";
      cerr << "Changed to " << _rows - 1 << endl;
      _rangeErrorCount--;
    }

    index = _rows - 1;
  }

  return _el[index];
}

//
//-------------------------// 
//
template <class Type>
void
Mat<Type>::operator () (const Mat<Type>& A)
{
  if ((A._maxrows != _maxrows) || (A._maxcols != _maxcols)) {
    _maxrows = A._maxrows;
    _maxcols = A._maxcols;
    _allocateEl();
  }

  _rows = A._rows;
  _cols = A._cols;

  if (_maxrows && _maxcols && _el)
    memcpy(*_el, *A._el, _maxrows*_maxcols*sizeof(Type));
}

//
//-------------------------// 
//


/**********************************************************************
Assignment operators for two dimensional matrix objects
***********************************************************************/

template <class Type>
Mat<Type>& 
Mat<Type>::operator = (const Mat<Type>& A)
{
  if (this == &A) 
    return *this;

  if ((A._maxrows != _maxrows) || (A._maxcols != _maxcols)) {
    _maxrows = A._maxrows;
    _maxcols = A._maxcols;
    _allocateEl();
  }
  
  _rows = A._rows;
  _cols = A._cols;

  if (_maxrows && _maxcols && _el)
    memcpy(*_el,*A._el, _maxrows*_maxcols*sizeof(Type));

  return *this;
}

template <class Type>
Mat<Type>& 
Mat<Type>::absorb(Mat<Type>& A)
{
  if (this == &A)
    return *this;

  // Delete current contents;
  if (_el) {
    if( _el[0] ) {
      delete [] _el[0];  // delete data block
    }
    delete [] _el; // delete row of pointers  
  }
  
  // Copy all from A
  _maxrows = A._maxrows;
  _maxcols = A._maxcols;
  _rows    = A._rows;
  _cols    = A._cols;
  _el      = A._el;
  
  // Empty A
  A._maxrows = A._maxcols = A._rows = A._cols = 0;
  A._el = 0;

  return *this;
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>&
Mat<Type>::pad(unsigned nrows, unsigned ncols, int row, int col, Type value)
{
  if ((nrows == _rows) && (ncols == _cols) && !row && !col)
    return *this;

  char tempFile[256];
  get_temp_filename(tempFile);

  // Save the current matrix to disk
  if (flushToDisk && saveRaw(tempFile)) {
    unsigned argRows = _rows;
    unsigned argCols = _cols;

    clear();

    _rows = _maxrows = nrows;
    _cols = _maxcols = ncols;

    _allocateEl();
    fill(value);

    insert(tempFile, argRows, argCols, row, col);
  }
  else {
    // Saving fails, create a new matrix
    Mat<Type> result(nrows, ncols, value);
    result.insert(*this, row, col);
    
    this->absorb(result);
  }

  unlink(tempFile);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::resize(unsigned nrows, unsigned ncols, int row, int col)
{
  if ((nrows == _rows) && (ncols == _cols))
    return *this;

  // Origin maintained, shrinking: no reallocation necessary
  if (!row && !col && (nrows <= _maxrows) && (ncols <= _maxcols)) {
    _rows = nrows;
    _cols = ncols;

    cerr << "This type of resizing is insecure!! Should be fixed..." << endl;

    return *this;
  }
  
  // This could be done by flushing to disk in order to save memory
  Mat<Type> T(nrows, ncols);

  Type *elPtr = T._el[0];
  for (int i = 0, i0 = row; i < nrows; i++, i0++) {
    Boolean valid = (i0 >= 0) && (i0 < _rows);
    for (int j = 0, j0 = col; j < ncols; j++, j0++)
      *elPtr++ = (valid && (j0 >= 0) && (j0 < _cols)) ? _el[i0][j0] : Type(0);
  }

  absorb(T);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::resize()
{
  _rows = _maxrows;
  _cols = _maxcols;
  
  cerr << "This type of resizing is insecure!! Should be fixed..." << endl;

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::appendRight(const Mat<Type>& A) const
{
  using std::max;		// (bert) force it to find the right max().

  Mat<Type> Temp(max(_rows, A._rows), _cols + A._cols);
  Temp.insert(*this);
  Temp.insert(A, 0, _cols);

  return Temp;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::appendBelow(const Mat<Type>& A) const
{
  using std::max;		// (bert) force it to use the right max()

  Mat<Type> Temp(_rows + A._rows, max(_cols, A._cols));
  Temp.insert(*this);
  Temp.insert(A, _rows, 0);

  return Temp;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::swapRows(unsigned r1, unsigned r2)
{
  if (r1 == r2)
    return *this;

  if ((r1 >= _rows) || (r2 >= _rows)) {
    cerr << "Error in swapRows: improper row indices " << r1 
	 << "," << r2 << " for matrix with " << _rows << " rows" << endl;

    return *this;
  }

  Type *r1ptr = _el[r1];
  Type *r2ptr = _el[r2];

  for (unsigned i = _cols; i; i--, r1ptr++, r2ptr++)
    swap(*r1ptr, *r2ptr);

  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::swapCols(unsigned c1, unsigned c2)
{
  if (c1 == c2)
    return *this;

  if ((c1 >= _cols) || (c2 >= _cols)) {
    cerr << "Error in swapCols: improper column indices " << c1 
	 << "," << c2 << " for matrix with " << _cols << " cols" << endl;

    return *this;
  }

  for (unsigned r = 0; r < _rows; r++)
    swap(_el[r][c1], _el[r][c2]);
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>
Mat<Type>::residual(unsigned row, unsigned col) const
{
  if ((_rows <= 1) || (_cols <= 1) || (row >= _rows) || (col >= _cols)) {
    cerr << "Error: residual(" << row << ", " << col << ") of " 
	 << _rows << "x" << _cols << " matrix." << endl;
    return *this;
  }

  Mat<Type> result(_rows - 1, _cols - 1);

  Type **rowPtr       = _el;
  Type **resultRowPtr = result._el;

  for (unsigned i = 0; i < _rows; i++) {
    Type *elPtr = *rowPtr++;
    if (i != row) {
      Type *resultElPtr = *resultRowPtr++;
      for (unsigned j = 0; j < _cols; j++) {
	if (j != col)
	  *resultElPtr++ = *elPtr;
	elPtr++;
      }
    }
  }
   
  return result;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::insert(const Mat<Type>& A, int row, int col)
{
  // This could be written more CPU-efficient
  const Type **ArowPtr = (const Type **) A._el;
  int          destRow = row;
  for (unsigned i = 0; i < A._rows; i++, destRow++) {
    Boolean valid = (destRow >= 0) && (destRow < _rows);
    const Type *AelPtr  = *ArowPtr++;
    int         destCol = col;
    for (unsigned j = 0; j < A._cols; j++, destCol++, AelPtr++)
      if (valid && (destCol >= 0) && (destCol < _cols))
	_el[destRow][destCol] = *AelPtr;
  }

  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::insert(const char *path, unsigned nrows, unsigned ncols, int row, int col)
{
  InputFile argFile(path);
  if (!argFile) {
    cerr << "Couldn't open file " << path << endl;
    return *this;
  }

  _checkMatrixDimensions(path, nrows, ncols);

  // Create a row-sized buffer
  Type *buffer = 0;
  allocateArray(ncols, buffer);
  if (!buffer) {
    cerr << "Couldn't allocate buffer" << endl;
    return *this;
  }

  int destRow = row;
  for (unsigned i = 0; i < nrows; i++, destRow++) {
    if (!(argFile.stream().read((char *) buffer, ncols*sizeof(Type)))) {
      cerr << "Error while reading file " << path << endl;
      freeArray(buffer);
      return *this;
    }
    Boolean valid = (destRow >= 0) && (destRow < _rows);
    const Type *bufferPtr = buffer;
    int         destCol = col;
    for (unsigned j = 0; j < ncols; j++, destCol++, bufferPtr++) {
      if (valid && (destCol >= 0) && (destCol < _cols))
	_el[destRow][destCol] = *bufferPtr;
    }
  }

  freeArray(buffer);
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fill(Type value)
{
  Type *elPtr = _el[0];
  for (unsigned i = _rows ; i ; i--)
    for (unsigned j = _cols ; j ; j--) 
      *elPtr++ = value ;
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fill(Type value, unsigned r1, unsigned r2, unsigned c1, unsigned c2)
{
  if ((r2 < r1) || (r2 >= _rows) || (c2 < c1) || (c2 >= _cols)) {
    cerr << "Error in Mat::fill: invalid row or column arguments." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }
    
  Type **rowPtr = _el + r1;
  for (unsigned i = r2 - r1 + 1 ; i ; i--) {
    Type *elPtr = *rowPtr++ + c1;  
    for (unsigned j = c2 - c1 + 1 ; j ; j--) 
      *elPtr++ = value ;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fillEllips(Type value, double rowDiameter, double colDiameter)
{
  double row = double(_rows - 1)/2;
  double col = double(_cols - 1)/2;

  if (rowDiameter <= 0)
    rowDiameter = _rows;
  if (colDiameter <= 0)
    colDiameter = _cols;

  double rowRadiusSqr = SQR(rowDiameter/2);
  double colRadiusSqr = SQR(colDiameter/2);

  Type *elPtr = _el[0];
  for (unsigned i = 0; i < _rows; i++) {
    double rowDist = SQR(i - row)/rowRadiusSqr;
    for (unsigned j = 0; j < _cols; j++, elPtr++) 
      if (rowDist + SQR(j - col)/colRadiusSqr <= 1)
	*elPtr = value;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::fillEllips(Type value, double row, double col, double rowDiameter, 
		      double colDiameter)
{
  if (rowDiameter <= 0)
    rowDiameter = 2*MIN(row + 0.5, _rows - row - 0.5);
  if (colDiameter <= 0)
    colDiameter = 2*MIN(col + 0.5, _cols - col - 0.5);

  double rowRadiusSqr = SQR(rowDiameter/2);
  double colRadiusSqr = SQR(colDiameter/2);

  Type *elPtr = _el[0];
  for (unsigned i = 0; i < _rows; i++) {
    double rowDist = SQR(i - row)/rowRadiusSqr;
    for (unsigned j = 0; j < _cols; j++, elPtr++) 
      if (rowDist +  SQR(j - col)/colRadiusSqr <= 1)
	*elPtr = value;
  }
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::eye()
{
  unsigned mindim = (_rows > _cols) ? _cols : _rows;
  
  Type **dataPtr = _el;
  
  for(unsigned j=0; j< _rows;j++)
    for(unsigned k=0;k<_cols;k++)
      dataPtr[j][k]=(Type) 0;
  
  for (unsigned i=0 ; i < mindim ; i++)
    dataPtr[i][i] = (Type)1.0;
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::randuniform(double min, double max)
{
  double range = max - min;

  Type *elPtr = _el[0];
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--)
      *elPtr++ = Type(drand48() * range + min);
  
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::randnormal(double mean, double std)
{
   Type *elPtr = _el[0];
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols ; j; j--)
       *elPtr++ = Type(gauss(mean, std));
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::hamming()
{
   //note len = _rows
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len - 1);
   double alpha = 0;
   
   for(unsigned i=0 ; i < len ; i++){
      _el[i][0] = (Type)(0.54 - 0.46*::cos(alpha));
      alpha += alphaIncrement;
   }
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::hanning()
{
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len-1);
   double alpha = 0;
   
   for(unsigned i=0 ; i < len; i++) {
      _el[i][0] = (Type)(0.5 - 0.5*::cos(alpha));
      alpha += alphaIncrement;
   }
   
   return *this;
}
//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::blackman()
{
   unsigned len=_rows;
   double alphaIncrement = 8.0*atan(1.0)/(len - 1);
   double alpha = 0.0;
   
   for(unsigned i=0 ; i < len ; i++) {
      _el[i][0] = (Type)(0.42 - 0.5*::cos(alpha) + 0.08*::cos((2*alpha)));
      alpha += alphaIncrement;
   }
   
   return *this;
}


//This function creates a diagnonal matrix from the existing row or column 
//vector object

template <class Type> 
Mat<Type> 
Mat<Type>::diag() const
{ 
//Ensure that calling object is a vector
  unsigned dimensions=0;

  if(this->isvector())
   { 
//Size of new square matrix will equal the size of the vector 
    dimensions=this->length();  
   }  
  else
    {
      cerr<<"Error:calling object is not a row or column vector"<<endl;
      exit(1);
    }
   
  Mat<Type> Diagonalized_Matrix(dimensions,dimensions,(Type)0);
  for(unsigned i=0 ; i < dimensions ; i++)
//The () operator for the calling object is used to access each element in order    
    Diagonalized_Matrix(i,i)=this->operator() (i);
   
  return(Diagonalized_Matrix);
}  

//This function creates a diagnonal matrix from the existing row or column 
//vector object

/****************************************************************************
This function creates a toeplitz matrix from the two input vectors

   ie.  if c=[a b c] and r=[A B C], the result would be  a B C
                                                         b a B
                                                         c b a
*****************************************************************************/
 
template<class Type>
Mat<Type>
Mat<Type>::toeplitz(const Mat<Type>& r) const
{
//Make sure that both input arguments are vectors
  if(!(r.isvector() && this->isvector() ))
   {
     cerr<<"Error:One or both of the input arguments is/are not a vector"<<endl;
     exit(1);
   }
//temporarily store all the lengths so we can safely and quickly use them in
//math
  unsigned c_length= this->length();
  unsigned r_length=r.length();
  unsigned Row_Buffer_length=(r_length - 1 + c_length);
  Mat<Type> Row_Buffer(1,Row_Buffer_length,(Type)0);

//Fill the row buffer with the two vectors concatenated
//The vector (calling object) c is concatenated backwards to the front of 
//the row buffer  

  for(unsigned i=0;i < c_length;i++)
    Row_Buffer(i)= this->operator () (c_length - 1 - i);
  for(unsigned j=1; j < r_length; j++)
    Row_Buffer(c_length - 1 + j)=r(j);
  
  Mat<Type> Toeplitz_M(1,r_length,(Type)0);
//Get first row of toeplitz matrix  
  Toeplitz_M=Row_Buffer(0,0,c_length -1,Row_Buffer_length -1);
//Construct the matrix by rolling and appending sections of Row_Buffer
  
  for(unsigned k=2;k <= c_length;k++)
    Toeplitz_M = Toeplitz_M.appendBelow(Row_Buffer(0,0,c_length-k,Row_Buffer_length-k));
  
  return(Toeplitz_M);   
} 

/*********************************************************************/
//Display functions for two dimensional matrix objects
/*********************************************************************/

/****************************************************/
//This function displays the entire matrice
//This one does not work for <type> complex elements
/***************************************************/

template <class Type>
ostream&
Mat<Type>::display(ostream& os) const
{
  for (unsigned i=0 ; i < _rows ; i++) {
    for (unsigned j=0 ; j < _cols ; j++)
      os << _el[i][j] << " ";
    os << endl;
  }
  return os;
}

/**********************************************************************
This function will display a subsection of the matrix determined by 
input arguments. The first argument is the starting row(the first row
is row 0) followed by the ending row.  The third and fourth argument
are analagous for the start and end columns. A seperate function is
used to display matrices with complex elements
***********************************************************************/

template <class Type>
ostream&
Mat<Type>::display(ostream& os, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const
{
//Make sure that the stop row and stop column are greater than the 
//start row and start column

  if ((r1 > r2) || (c1 > c2)){
    cerr << "Error in display: improper row or column sizes." << endl;
    cerr << r1 << " to " << r2 << " and" << endl;
    cerr << c1 << " to " << c2 << endl;
    exit(1);
  }

//Check to see if the stop row and stop column are exclusively less than
//the actual number of rows and columns in the matrix.  This
//prevents accidental access to an undefined memory location

  if ((r2 >=_rows) || (c2 >=_cols)) {
    cerr<<"The requested _rows or columns are not defined for this ";
    cerr<<"matrix"<<endl;
    exit(1);
  } 
  
  for (unsigned i=r1 ; i <= r2 ; i++) {
    for (unsigned j=c1 ; j <= c2 ; j++)
      os << _el[i][j] << " ";
    os << endl;
  }
  return os;
}

/**********************************************************************
This function will display a subsection of the matrix determined by 
input arguments.  The first argument is the starting row(the first row
is row 0) followed by the ending row.  The third and fourth argument
are analagous for the start and end columns.  This works or only matrices
with complex elements.
***********************************************************************/

//
//------------
//
template <class Type>
int
Mat<Type>::operator != (const Mat<Type>& Arg) const
{
  if ((_rows != Arg._rows) || (_cols != Arg._cols))
    return 1;

  Type **argRowPtr = Arg._el;
  Type **rowPtr    = _el;
  for (unsigned i=_rows ; i != 0 ; i--) {
    Type *argColPtr = *argRowPtr++;
    Type *colPtr    = *rowPtr++;
    for (unsigned j=_cols ; j != 0 ; j--)
      if (*colPtr++ != *argColPtr++)
	return 1;
  }
  
  return 0;
}

template <class Type>
Mat<Type>& 
Mat<Type>::operator += (dcomplex x)
{
  Type addend = Type(real(x));

  Type *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr += addend;
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::operator *= (dcomplex x)
{
  double scale = real(x);

  Type *elPtr = _el[0];
  for (unsigned i=_rows ; i != 0 ; i--)
    for (unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr = Type(scale * *elPtr);
  
  return *this;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::t() const            // Transpose operation
{
   Mat<Type> Temp(_cols, _rows);
   Type **tempDataRowPtr = Temp._el;
   for (unsigned i=0 ; i < _cols ; i++) {
      Type *tempDataColPtr = *tempDataRowPtr++;
      for (unsigned j=0 ; j < _rows ; j++)
	 *tempDataColPtr++ = _el[j][i];
   }
   
   return Temp;
}

//************************
//reverse the filter or rotate by 180
//used for convolution and morphology 
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::rotate180() const           
{
   Mat<Type> Temp(_rows,_cols);
   unsigned kr,kc;   
   for (unsigned i=0 ; i < _rows ; i++) {
      kr = _rows -i -1;
      for (unsigned j=0 ; j < _cols ; j++){
	 kc = _cols -j -1;
	 Temp(kr,kc) = _el[i][j];
      }
   }
   
   return Temp;
}

/**************************

inv - matrix inversion - Gauss-Jordan elimination with full pivoting.

The returned matrix is always a Type matrix.  Exits and prints
error message with singular matrices or bad MATRIX structures.

***************************/
template <class Type>
Mat<Type> 
Mat<Type>::inv() const
{
//Good results if used on float or double only
// check for square matrix 
  if(_rows != _cols) {
    cerr << endl << "Mat inversion ERROR: non-square, size = " 
	 << _rows << " x " << _cols << endl;
    return Mat<Type>();
  }
   
// check pointer 
  if(!_el) {
    cerr << endl << "ERROR: No matrix for matrix invert" << endl;
    return Mat<Type>();
  }
  
// size of square matrix 
  unsigned      n = _rows;
  Mat<Type>   Ai(*this);
  Type **a = Ai._el;
  
// allocate index arrays and set to zero 
  unsigned *pivot_flag = (unsigned *) calloc(n,sizeof(unsigned));
  unsigned *swap_row = (unsigned *) calloc(n,sizeof(unsigned));
  unsigned *swap_col = (unsigned *) calloc(n,sizeof(unsigned));
  
  if(!pivot_flag || !swap_row || !swap_col) {
    cerr << "Allocation error in matrix invert" << endl;
    return Mat<Type>();
  }
  
  for(unsigned i = 0 ; i < n ; i++) {   // n iterations of pivoting 
      
// find the biggest pivot element 
    unsigned irow,icol;
    unsigned row, col;
    Type big = (Type)0;
    for(row = 0 ; row < n ; row++) {
      if(!pivot_flag[row]) {       // only unused pivots 
	for(col = 0 ; col < n ; col++) {
	  if(!pivot_flag[col]) {
	    Type abs_element = (Type)fabs((double)a[row][col]);
	    if(abs_element >= big) {
	      big = abs_element;
	      irow = row;
	      icol = col;
	    }
	  }
	}
      }
    }
    pivot_flag[icol]++;    // mark this pivot as used 
    
    // swap rows to make this diagonal the biggest absolute pivot 
    Ai.swapRows(irow, icol);

    /* Old code for row swap
    if(irow != icol) {
      for(unsigned col = 0 ; col < n ; col++) {
	Type temp = a[irow][col];
	a[irow][col] = a[icol][col];
	a[icol][col] = temp;
      }
    }
    */
    
    // store what we swaped 
    swap_row[i] = irow;
    swap_col[i] = icol;
      
    // bad news if the pivot is zero 
    if(a[icol][icol] == (Type) 0) {
      cerr << "ERROR: matrix invert: SINGULAR MATRIX" << endl;
      return Mat<Type>();
    }
    
    // divide the row by the pivot 
    //note here if the opreration is integer results are bad
    Type pivot_inverse = (Type)1/a[icol][icol];
    a[icol][icol] = (Type) 1;        // pivot = 1 to avoid round off 
    for(col = 0 ; col < n ; col++)
      a[icol][col] = a[icol][col]*pivot_inverse;
    
    // fix the other _rows by subtracting 
    for(row = 0 ; row < n ; row++) {
      if(row != icol) {
	Type temp = a[row][icol];
	a[row][icol] = (Type) 0;
	for( col = 0 ; col < n ; col++)
	  a[row][col] = a[row][col]-a[icol][col]*temp;
      }
    }
  }
  
  // fix the affect of all the swaps for final answer 
  for(int swap = n-1 ; swap >= 0 ; swap--) {
    Ai.swapCols(swap_row[swap], swap_col[swap]);

    /* Old code for column swap
    if(swap_row[swap] != swap_col[swap]) {
      for(unsigned row = 0 ; row < n ; row++) {
	Type temp = a[row][swap_row[swap]];
	a[row][swap_row[swap]] = a[row][swap_col[swap]];
	a[row][swap_col[swap]] = temp;
      }
    }
    */
  }

     
  // free up all the index arrays 
  free((char *)pivot_flag);
  free((char *)swap_row);
  free((char *)swap_col);
  
  return Ai;
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::h() const        // Hermitian Transpose operation
{
   Mat<Type> Temp(_cols, _rows);
   Type **tempDataRowPtr = Temp._el;
   for (unsigned i=0 ; i < _cols ; i++) {
      Type *tempDataColPtr = *tempDataRowPtr++;
      for (unsigned j=0 ; j < _rows ; j++)
	 *tempDataColPtr++ = _el[j][i];
   }
   return Temp;
}

//
//-------------------------// 
//
template <class Type>
Type 
Mat<Type>::min(unsigned *row, unsigned *col) const
{
   Type min = **_el;
   unsigned    rowOfMin = 0, colOfMin = 0;
   unsigned    rowindx = 0, colindx = 0;

   Type **rowPtr = _el;
   for(unsigned i = _rows ; i != 0 ; i--) {
      Type *colPtr = *rowPtr++;
      colindx = 0;
      for(unsigned j = _cols ; j != 0 ; j--){
	 if (*colPtr < min) {
	    min = *colPtr;
	    rowOfMin = rowindx;
	    colOfMin = colindx;
	 }
	 colPtr++;
         colindx++; 
      }
      rowindx++;
   }
   
   if (row)
      *row = rowOfMin;
   if (col)
      *col = colOfMin;
   
   return(min);
}

//
//-------------------------// 
//
template <class Type>
Type 
Mat<Type>::max(unsigned *row, unsigned *col) const
{
   Type max = **_el;
   unsigned    rowOfMax = 0, colOfMax = 0;
   unsigned    rowindx = 0, colindx = 0;
   
   Type **rowPtr = _el;
   for(unsigned i = _rows ; i != 0 ; i--) {
      Type *colPtr = *rowPtr++;
      colindx = 0;
      for(unsigned j = _cols ; j != 0 ; j--){
	 if (*colPtr > max) {
	    max = *colPtr;
	    rowOfMax = rowindx;
	    colOfMax = colindx;
	 }
	 colPtr++;
         colindx++;
      }
      rowindx++;
   }
   
   if (row)
      *row = rowOfMax;
   if (col)
      *col = colOfMax;
   
   return(max);
}

//
//-------------------------// 
//
template <class Type>
Type
Mat<Type>::median(Type minVal, Type maxVal) const
{
  return array(minVal, maxVal).medianVolatile();
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::csum() const
{
  dcomplex result = 0;
   
  Type **rowPtr = _el;
  for(unsigned i = _rows ; i != 0 ; i--) {
    Type *colPtr = *rowPtr++;
    for(unsigned j = _cols ; j != 0 ; j--)
      result += *colPtr++;
  }
  
  return result;
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::csum2() const
{
  dcomplex result = 0;
   
  Type **rowPtr = _el;
  for(unsigned i = _rows ; i != 0 ; i--) {
    Type *colPtr = *rowPtr++;
    for(unsigned j = _cols ; j != 0 ; j--, colPtr++)
      result += dcomplex(*colPtr * *colPtr);
  }
  
  return result;
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::ctrace() const
{
   unsigned mindim = (_rows >= _cols) ? _cols : _rows;
   
   Type **dataPtr = _el;
   dcomplex temp = 0;
   for (unsigned i=0 ; i < mindim ; i++)
      temp += dataPtr[i][i];
   
   return(temp);
}

//
//-------------------------// 
//
template <class Type>
dcomplex
Mat<Type>::cdet() const
{
  dcomplex determinant = 0;

  if (!_rows || (_rows != _cols)) {
    cerr << "Error: determinant of non-square or empty matrix" << endl;
    return determinant;
  }

  if (_cols > 1) {
    int factor = 1;
    Type *elPtr = *_el;
    for(unsigned col = 0; col < _cols; col++, factor *= -1)
      determinant += dcomplex(*elPtr++) * (residual(0, col).det() * factor);
  }
  else 
    determinant = **_el;

  return determinant; 
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyElementWise(double (*function)(double))
{
  Type *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (Type) function((double) *elPtr);
  
  return(*this);
}

template <class Type>
Mat<Type>&
Mat<Type>::applyElementWiseC2D(double (*function)(dcomplex))
{
   Type *elPtr = _el[0];
   
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols; j; j--, elPtr++)
       *elPtr = (Type) function((dcomplex) *elPtr);
   
   return(*this);
}

template <class Type>
Mat<Type>&
Mat<Type>::applyElementWiseC2C(dcomplex (*function)(dcomplex))
{
   Type *elPtr = _el[0];
   
   for(unsigned i = _rows; i; i--)
     for(unsigned j = _cols; j; j--, elPtr++)
       *elPtr = Type(real(function((dcomplex) *elPtr)));
   
   return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyIndexFunction(IndexFunction F)
{
  Type *elPtr = *_el;
  for (unsigned r = 0; r < _rows; r++)
    for (unsigned c = 0; c < _rows; c++)
      *elPtr++ = Type(F(r, c));

  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::applyIndexFunction(ComplexIndexFunction F)
{
  Type *elPtr = *_el;
  for (unsigned r = 0; r < _rows; r++)
    for (unsigned c = 0; c < _rows; c++)
      *elPtr++ = Type(asDouble(F(r, c)));

  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::exp()
{
  return applyElementWise(::exp);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::log()
{
  return applyElementWise(::log);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::cos()
{
  return applyElementWise(::cos);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::sin()
{
  Type *elPtr = _el[0];
  
  for(unsigned i = _rows; i; i--)
    for(unsigned j = _cols; j; j--, elPtr++)
      *elPtr = (Type) ::sin(double(*elPtr));
  
  return(*this);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::round()
{
  return applyElementWise(::rint);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type>&
Mat<Type>::abs()
{
  return applyElementWise(::fabs);
}

//
//-------------------------// 
//

template <class Type>
Mat<Type>&
Mat<Type>::sqrt()
{
  return applyElementWise(::sqrt);
}

template <class Type>
Mat<Type>&
Mat<Type>::pow(double exponent)
{
  Type *elPtr = _el[0];
  for(unsigned i=_rows ; i != 0 ; i--)
    for(unsigned j=_cols ; j != 0 ; j--, elPtr++)
      *elPtr = (Type)(::pow((double)*elPtr, exponent));
  
  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::conj()
{
  return *this;
}

//
//-------------------------//
//
template <class Type>
Mat<Type>
Mat<Type>::transposeXself() const
{
  Mat<Type>  result(_cols, _cols);
  Type     **resultEl = result._el;

  for (unsigned i = 0; i < _cols; i++) {
    Type *resultPtr = resultEl[i];
    
    for (unsigned j = 0; j < i; j++) {
      Type& value = *resultPtr++;
      value = Type(0);
      for (unsigned k = 0; k < _rows; k++)
	value += _el[k][i] * _el[k][j];

      resultEl[j][i] = value;
    }
    
    Type& value = *resultPtr++;
    value = Type(0);
    for (unsigned k = 0; k < _rows; k++)
      value += _el[k][i] * _el[k][i];
  }
  
  return result;
}

//
//-------------------------//
//
template <class Type>
void 
Mat<Type>::eig(Mat<Type>& D, Mat<Type>& V) const
{
   //Mat<double> Ad = (DblMat) (this);
   
   double  *d;
   unsigned np1;
   unsigned *positions, itemp;
   double  temp;
   
   // if (!Ad.getEl()) {
   if (!_el) {
      printf("eig: invalid input matrix pointer\n");
      exit(1);
   }
   
   if (_rows != _cols) {
      cerr << "eig: matrix is not square -- " << _rows << " x " << _cols << endl;
      exit(1);
   }
   
   np1 = _rows + 1;
   Mat<double> A(np1,np1);
   Mat<double> v(np1,np1);
   d = (double *) malloc((np1)*sizeof(double));
   unsigned i, j;
   for (i=1 ; i <= _rows ; i++)
      for (j=1 ; j <= _cols ; j++)
	 A(i,j) = (double)_el[i-1][j-1];
   //A(i,j) = Ad.getEl[i-1][j-1];
   
   jacobi((double **) A.getEl(), _rows, d,(double **) v.getEl());
   
// sort eigenvalues in descending order 
   
   positions = (unsigned *) malloc( (_rows + 1) * sizeof(unsigned));
   for (i=1 ; i <= _rows ; i++)
      positions[i] = i;
   for (i=1 ; i <= _rows ; i++)
      for (j=1 ; j < _rows ; j++){
	 if (d[j] < d[j+1]) {
	    temp = d[j]; d[j] = d[j+1]; d[j+1] = temp;
	    itemp = positions[j]; 
	    positions[j] = positions[j+1]; // Keep track of original positions 
	    positions[j+1] = itemp;
	 } 
      }
   
   for (i=0 ; i < _rows ; i++)
      for (j=0 ; j < _cols ; j++) {
	 V(i,j) = (Type)v(i+1,positions[j+1]);
	 D(i,j) = (Type)0.0;
      }
   for (i=0 ; i < _rows ; i++)
      D(i,i) = (Type)d[i+1];
   
   
   free((char *)d);
   free((char *)positions);
}

//
//-------------------------// 
//
template <class Type>
Mat<Type> 
Mat<Type>::house() const
{
   if (_cols != 1) {
      cerr << "Error: input to house is not a column vector." << endl;
      exit(1);
   }
   
   unsigned n = _rows;
   Type mu = (Type) 0.0;
   Type **rowPtr = _el;
   for(unsigned i=n ; i != 0 ; i--, rowPtr++) {
      mu += (**rowPtr) * (**rowPtr);
   }
   
   mu = (Type)(::sqrt((double)mu));
   Mat<Type> v(*this);
   if (mu != (Type) 0.0){
      Type xFirst = **_el;
      Type beta = (xFirst > (Type) 0.0) ? xFirst + mu : xFirst - mu;
      
      Type **rowPtr  = _el + 1;
      Type **vRowPtr = v._el + 1;
      for(unsigned i=n-1 ; i != 0 ; i--, rowPtr++, vRowPtr++){
	 **vRowPtr = **rowPtr / beta;
      }
   }
   
   **v._el = (Type) 1.0;
   
   return(v);
}

//
//-------------------------//Not sure 
//
template <class Type>
Mat<Type> 
Mat<Type>::rowhouse(const Mat<Type>& V) const
{
   Type beta;

   if (V._cols != 1) {
      cerr << "Error: input to rowhouse is not a column vector." << endl;
      exit(1);
   }
   
   if (V._rows != _rows) {
      cerr << "Error: vector input to rowhouse is wrong length." << endl;
      exit(1);
   }
   
   Type nv = (Type) 0.0;
   Type **vRowPtr = V._el;
   for (unsigned i=V._rows ; i != 0 ; i--, vRowPtr++)
      nv += (**vRowPtr) * (**vRowPtr);
   
   if (nv == (Type) 0.0){
      cerr << "Error: vector input to rowhouse is all Zeros." << endl;
      exit(1);
   }
   beta = (Type)(-2.0)/nv;
   Mat <Type> w(_cols,1);
   w =    (*this).t();
   w =    w * V;
   w =  w * beta;
   return ( *this + V * w.t() );
}


//
//-------------------------//
//
template <class Type>
Mat<Type>
Mat<Type>::convolv2d(const Mat<Type>& filter) const
{
  unsigned    row,column, i, j;
  unsigned    dead_rows, dead_cols, kernel_rows, kernel_cols, l_cols, l_rows;
  Type   *outptr, *inptr, *filterptr;

  Mat<Type> out(_rows,_cols);
  Mat<Type> filterrev = filter.rotate180();
  kernel_rows = filter.getrows();
  kernel_cols = filter.getcols();
  dead_rows = kernel_rows/2;
  dead_cols = kernel_cols/2;
  l_rows = 2*dead_rows + _rows;
  l_cols = 2*dead_cols +_cols;
  Mat<Type> large(l_rows, l_cols);
  Mat<Type> outlarge(l_rows, l_cols);

/*copy in into center of large*/
  for(i=0; i < _rows ; i++){
    for(j=0; j < _cols ; j++)
      large._el[i+dead_rows][j+dead_cols] = _el[i][j];
  }

  for(row=0; row < large._rows -kernel_rows+1; row++){
    outptr = outlarge._el[row+dead_rows];
    for (column=0; column<large._cols-kernel_cols+1; column++){
      Type sumval = 0;
      for (i=0; i< kernel_rows; i++){
	inptr = large._el[i+row];
	inptr = inptr + column;
	filterptr = (Type *) filterrev._el[i];
	for (j=0; j< kernel_cols; j++)
	  sumval += (*inptr++) * (*filterptr++);
      }

      outptr[column+dead_cols] = sumval;
    }
  }
  
  for(i=0; i < _rows ; i++){
    for(j=0; j < _cols ; j++)
      out._el[i][j] = outlarge._el[i+dead_rows][j+dead_cols];
  }

  return(out);
}

/***************/

template <class Type>
Histogram
Mat<Type>::histogram(double minin, double maxin, unsigned n) const
{
  if (maxin <= minin) {
    minin = min();
    maxin = max();
  }

  Histogram histogram(minin, maxin, n);

  const Type **rowPtr = (const Type **) _el;
  for (unsigned i = _rows; i; i--) {
    const Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--)
      histogram.add(*elPtr++);
  }
  
  return histogram;
}

/*
//this histogram works in the following way
//if there is no input in the histogram ex: histogram(), then
//it takes the defaults 0,0 and does the histogram from the smallest
//value in the image to the largest value.
//if the input is histogram(0,255) it will give the histogram
//between 0 and 255.  If the input is 10,20 it will give only 10,20
//values

template <class Type>
Mat<Type>
Mat<Type>::histogram(int minin, int maxin) const
{
  if(minin == maxin){
    minin=(int)min();
    maxin=(int)max();
  }
  int range = maxin - minin + 1;
  Mat<Type> histogram1;
  histogram1 = Zeros<Type>(1, range);
  Type *histPtr = histogram1._el[0];
  Type **rowPtr = _el;
  for (int i = _rows; i != 0; i--) {
    Type *colPtr = *rowPtr++;
    for (int j = _cols; j != 0; j--, colPtr++) {
      if ((*colPtr >= minin) && (*colPtr <= maxin)) {
	int roundedEl = (int) rint((double)*colPtr);
	(*(histPtr + roundedEl))++;
      }
    }
  }
  
  return(histogram1);
}
*/

//
//-------------------------//
//
template <class Type>
Mat<Type>&
Mat<Type>::histmod(const Histogram& hist1, const Histogram& hist2)
{
  return map(equalize(hist1, hist2));
}

//
//-------------------------//
//
template <class Type>
SimpleArray<Type>
Mat<Type>::array(Type minVal, Type maxVal) const
{
  unsigned N = 0;
  Type   **rowPtr;
  Type    *elPtr;

  if (maxVal <= minVal) {
    minVal = min();
    maxVal = max();
    N = nElements();
  }
  else {
    // Scan matrix to determine # elements in the range
    rowPtr = _el;
    for (unsigned i = _rows; i; i--) {
      elPtr = *rowPtr++;
      for (unsigned j = _cols; j; j--, elPtr++)
	if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	  N++;
    }
  }
  
  SimpleArray<Type> array(N);

  if (N) {
    Type *arrayPtr = array.contents();
    rowPtr = _el;
    for (unsigned i = _rows; i; i--) {
      elPtr = *rowPtr++;
      for (unsigned j = _cols; j; j--)
	if ((*elPtr >= minVal) && (*elPtr <= maxVal))
	  *arrayPtr++ = *elPtr++;
    }
  }

  return array;
}

//
//-------------------------//
//
template <class Type>
void 
Mat<Type>::qr(Mat<Type>& R, Mat<Type>& Q) const
{
   /*****************Don't need this because R and Q get reassigned later
   if ((_rows != R._rows) && (_cols != R._cols)) {
      cerr << "Error: R matrix of incompatible size for QR." << endl;
      exit(1);
   }
   if ((_rows != Q._rows) && (_rows != Q._cols)) {
      cerr << "Error: Q matrix of incompatible size for QR." << endl;
      exit(1);
   }
   ***********************************************************************/
   
   Mat<Type> Atmp(*this);
   //Mat<Type> vtmp(Zeros<Type>(_rows,1));
   //This vector vtmp should automatically be initialized to zero
   Mat<Type> vtmp(_rows,1);
   Mat<Type> QTemp(_rows,_rows);
//This below instantiation used to be at position 12345
   Mat<Type> Temp(*this);
   unsigned mindim = (_rows > _cols) ? _cols : _rows;
   
// This can be made more efficient -- Alex, 2/4/94

   unsigned i, k;
   int      j; 
   for(j=0 ; j < mindim ; j++){
      vtmp.resize(_rows-j,1);
      
      for(i=j ; i < _rows ; i++)
	 vtmp(i-j,0) = Atmp(i,j);
      
      vtmp = vtmp.house();
//Position 12345
      Temp(Temp.rowhouse(vtmp));
      
      for(i=j ; i < _rows ; i++)
	 for(k=j ; k < _cols ; k++)
	    Atmp(i,k) = Temp(i-j,k-j);
      
      if (j < (_rows-1)){
	 for(i=j+1 ; i < _rows ; i++)
	    Atmp(i,j) = vtmp(i-j,0);
      }
      
      Temp.resize(_rows-j-1,_cols-j-1);
      for(i=0 ; i < _rows-j-1 ; i++)
	 for(k=0 ; k < _cols-j-1 ; k++)
	    Temp(i,k) = Atmp(i+j+1,k+j+1);
   }
   
   R = Mat<Type>(_rows,_cols);
   for(i=0 ; i < mindim ; i++) 
      for(k=i ; k < mindim ; k++)
	 R(i,k) = Atmp(i,k);
   
   vtmp.resize();
   Q = Mat<Type>(_rows,_rows);
   Q.eye();
   for(j = mindim - 1 ; j >= 0 ; j--){
      vtmp.resize(_rows-j,1);
      QTemp.resize(_rows-j,_rows-j);
      vtmp(0,0) = (Type) 1;
      for(i=1 ; i < _rows-j ; i++)
	 vtmp(i,0) = Atmp(i+j,j);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    QTemp(i,k) = Q(i+j,k+j);
      QTemp = QTemp.rowhouse(vtmp);
      for( i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    Q(i+j,k+j) = QTemp(i,k);
   }
   vtmp.resize();
   if (_rows < _cols){
      for(j=_rows ; j < _cols ; j++){
	 for(i=0 ; i < _rows ; i++){
	    vtmp(i,0) = _el[i][j];
	 }
	 vtmp = Q.t() * vtmp;
	 for(i=0 ; i < _rows ; i++){
	    R(i,j) = vtmp(i,0);
	 }
      }
   }
}

/*********************************************************
template <class Type>
void 
Mat<Type>::qrtwo(Mat<Type>& R, Mat<Type>& Q) const
{
   cerr << "Dcomplex 3 qr not implemented" << endl;
}
**********************************************************/
//
//-------------------------//
//
/*****************************************************************************
This qr function only works with the _rows greater than or equal to the number 
of columns.  It may be possible to modify this so that it works with any size
matrix.
******************************************************************************/
template <class Type>
void 
Mat<Type>::qr(Mat<Type>& R, Mat<Type>& Q, Mat<Type>& P) const
{
// This can be made more efficient -- Alex, 2/4/94 
   
   int i, j, k, r, *piv, mindim;
   unsigned rloc, cloc;
   Type tau, dtmp;
   Mat<Type> Atmp(*this);

/*********changed 6/17/94******************   
   P = Eye<Type>(_cols);
   R = Zeros<Type>(_rows,_cols);
*******************************************/
   if(_cols>_rows){
      cerr<<"_Rows must be greater than or equal to columns"<<endl;
      exit(1);
   } 
   P=Mat<Type>(_cols,_cols);
   P.eye();
   R=Mat<Type>(_rows,_cols);

   piv = (int *)malloc(_cols * sizeof(int));
   if (!piv){
      cerr << "Error forming piv integer vector in qr." << endl;
      cerr << "_rows = " << _rows << " _cols = " << _cols << endl;
      exit(1);
   }
   Mat<Type> c(_cols,1);
   Mat<Type> v(_rows,1);
   Mat<Type> tmp(_cols,1);
   Mat<Type> Temp(_rows,_cols);
   Mat<Type> QTemp(_rows,_rows);
   Mat<Type> PTemp(_cols,_cols);
   for(j=0 ; j < _cols ; j++){
      c(j) = (Type)0.0;
      piv[j] = 0;
      for(i=0 ; i < _rows ; i++)
	 c(j) = c(j) + Atmp(i,j)*Atmp(i,j);
   }
   r = -1;
   tau = c.max(&rloc, &cloc);
   k = rloc;

   
   while(tau > 0){
      r++;
      piv[r] = k;
      //Swap the rth and kth columns of A      
      for(i=0 ; i < _rows ; i++){
	 dtmp=Atmp(i,r); Atmp(i,r)=Atmp(i,k); Atmp(i,k)=dtmp;
      } 	 
      //Swap c(r) and c(k)      
      dtmp=c(r);   c(r)=c(k);     c(k)=dtmp;

        
/**********changed 6/17/94******************      
      PTemp = Eye<Type>(_cols);
********************************************/
//Obtain an identity matrix to perform swaps to calulate the Pth permutation
//matrix

      PTemp=Mat<Type> (_cols,_cols);
      PTemp.eye();

//swap row j=r and piv[j]=k of PTemp to create the jth permutation matrix
      PTemp(r,r)=(Type)0.0; 
      PTemp(k,k)=(Type)0.0; 
      PTemp(r,k)=(Type)1.0;
      PTemp(k,r)=(Type)1.0;
//Accumulate the Permutation matrix initially, P is an identity matrix
      P = P*PTemp;
//Resize Temp and v to prepare for the rth householder       
      Temp.resize(_rows-r,1);  
      v.resize(_rows-r,1);
//Temp gets elements A(r:m,r) 
      for(i=r ; i < _rows ; i++)
	 Temp(i-r,0) = Atmp(i,r);
//v gets house(A(r:m,r))
      v = Temp.house();
//Resize Temp to allow the rowhouse 
      Temp.resize(_rows-r,_cols-r);
//Temp=A(r:m,r:n)
      for(i=r ; i < _rows ; i++)
	 for(j=r ; j < _cols ; j++)
	    Temp(i-r,j-r) = Atmp(i,j);
//Temp=rowhouse(A(r:m,r:n),v(r:m))
      Temp = Temp.rowhouse(v);
//A(r:m,r:n)=Temp
      for(i=r ; i < _rows ; i++)
	 for(j=r ; j < _cols ; j++)
	    Atmp(i,j) = Temp(i-r,j-r);
//A(r+1:m,r)=v(r+1:n)
      for(i=r+1 ; i < _rows ; i++)
	 Atmp(i,r) = v(i-r);

      for(i=r+1 ; i < _cols ; i++)
	 c(i) = c(i) - Atmp(r,i)*Atmp(r,i);

      if ( r < (_cols -1))
      {
         /******************changed 6/21/94*************************
	 tmp.resize(_cols-1-r,1);
	 for(i=r+1 ; i < _cols ; i++)
	    tmp(i-r-1,0) = c(i);
	 tau = tmp.max(&rloc, &cloc);
         ***********************************************************/
         //tau gets max(c(r+1:n))
         tau = (c(r+1,_cols-1,0,0)).max(&rloc,&cloc); 
	 k = rloc + r + 1;
      }
      else
	 tau = (Type) 0;
   }
//end while
//R gets upper triangular part of A
   for(i=0 ; i < _rows ; i++)
      for(j=i ; j < _cols ; j++)
	 R(i,j) = Atmp(i,j);
   
   mindim = r+1;
   v.resize();
  
  /**************changed 6/17/94******************   
   Q = Eye<Type>(_rows);
   ***********************************************/
   Q=Mat<Type>(_rows,_rows);
   Q.eye();

   for( j = mindim - 1 ; j >= 0 ; j--){
      v.resize(_rows-j,1);
      QTemp.resize(_rows-j,_rows-j);
      v(0,0) = (Type) 1;
      for(i=1 ; i < _rows-j ; i++)
	 v(i,0) = Atmp(i+j,j);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    QTemp(i,k) = Q(i+j,k+j);
      QTemp = QTemp.rowhouse(v);
      for(i=0 ; i < _rows-j ; i++)
	 for(k=0 ; k < _rows-j ; k++)
	    Q(i+j,k+j) = QTemp(i,k);
   }
   v.resize();
   if (_rows < _cols){
      for(j=_rows ; j < _cols ; j++){
	 for(i=0 ; i < _rows ; i++){
	    v(i,0) = _el[i][j];
	 }
	 v = Q.t() * v;
	 for(i=0 ; i < _rows ; i++){
	    R(i,j) = v(i,0);
	 }
      }
   }
   free((char *)piv);
}

//
//-------------------------//
//
/*****************************************************************************
This svd function only works with the _rows greater than or equal to the number 
of columns.  It may be possible to modify this so that it works with any size
matrix.
******************************************************************************/

template <class Type>
void 
Mat<Type>::svd(Mat<Type>& U, Mat<Type>& S, Mat<Type>& V) const
{
// This can be made more efficient -- Alex, 2/4/94
/************************changed 6/27/94*************************************   
   if ((U._rows != _rows) || (U._cols != _rows)){
      cerr << "Error in svd: left eigenvector input matrix of wrong size." << endl;
      cerr << _rows << "x" << _rows << endl;
      exit(1);
   }
   if ((S._rows != _rows) || (S._cols != _cols)){
      cerr << "Error in svd: singular value input matrix of wrong size." << endl;
      cerr << _rows << "x" << _rows << endl;
      exit(1);
   }
   if ((V._rows != _cols) || (V._cols != _cols)){
      cerr << "Error in svd: right eigenvector input matrix of wrong size." << endl;
      cerr << _rows << "x" << _cols << endl;
      exit(1);
   }
******************************************************************************/
   
   if(_cols>_rows){
     cerr<<"_Rows must be greater than or equal to columns"<<endl;
     exit(1);
   } 
   
   U=Mat<Type> (_rows,_rows);
   S=Mat<Type> (_rows,_cols);
   V=Mat<Type> (_cols,_cols);
   
   Mat<Type> C(_cols,_cols);
   Mat<Type> D(_cols,_cols);
   Mat<Type> R(_rows,_cols);
   Mat<Type> P(_cols,_cols);
   
   C = t() * (*this);
   C.eig(D, V);
   
   ((*this)*V).qr(R,U,P);
   
   unsigned mindim = (_rows > _cols) ? _cols : _rows;
   
   S = U.t()*(*this)*V*P;
   unsigned i, j;
   for(i=0 ; i < mindim ; i++){
     Type& sii = S(i, i);
     if (double(sii) < 0.0) {
       sii = int(fabs(double(sii)));
       for(j=0 ; j < _rows ; j++)
	 U(j,i) = -U(j,i);
     }
   }
   V = V*P;
   
   for(i=0 ; i < _rows ; i++)
     for(j=0 ; j < _cols ; j++){
       if ( i != j)
	 S(i,j) = (Type)0.0;
     }
}

//
//-------------------------//
//
template < class Type>
Mat<Type>&
Mat<Type>::clip(Type minVal, Type maxVal, Type minFill, Type maxFill)
{
  Type **rowPtr = _el;
  for(unsigned i = _rows; i; i--) {
    Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--, elPtr++) {
      if (*elPtr < minVal)
	*elPtr = minFill;
      if (*elPtr > maxVal)
	*elPtr = maxFill;
    }
  }
  
  return(*this);
}

//
//-------------------------//
//
template < class Type>
Mat<Type>&
Mat<Type>::map(const ValueMap& valueMap)
{
  Type **rowPtr = _el;
  for(unsigned i = _rows; i; i--) {
    Type *elPtr = *rowPtr++;
    for (unsigned j = _cols; j; j--, elPtr++)
      *elPtr = Type(valueMap(double(*elPtr)));
  }
  
  return(*this);
}

//
//-------------------------//
//
//The scale function will scale the image between 0 and 255
//unless given a range with a minimum and maximum value
//then it finds 256 levels and adds the minimum to all values
 
//note that: ((x-min)*(Max-Min)/(max-min))+Min
//          = x*((Max-Min)/(max-min)) + ((Minmax - Maxmin)/(max-min))

template < class Type>
Mat<Type>&
Mat<Type>::scale(double minout, double maxout, double minin, double maxin)
{
  if (minin >= maxin) {
    minin = asDouble(min());
    maxin = asDouble(max());
  }

  map(LinearMap(minin, maxin, minout, maxout));
  clip(Type(minout), Type(maxout), Type(minout), Type(maxout));

  return *this;
}

template <class Type>
void 
Mat<Type>::linearinterpolate(unsigned Urow, unsigned Drow, unsigned Ucol, unsigned Dcol, Mat<Type>& y) const
{

   unsigned     i, nrows, ncols, out_rows, out_cols;
   unsigned     n, m;
   
   nrows = _rows;
   ncols = _cols;
   
   out_cols = ((ncols-1)*Urow+1)/Drow;
   out_rows = ((nrows-1)*Ucol+1)/Dcol;

/**************changed 6/23/94**********************   
   remake(y, out_rows, out_cols);
   y.fill((Type)0);
   y = Zeros<Type>(out_rows, out_cols);
   
   yi.fill((Type)0);
   yi = Zeros<Type>(nrows, out_cols);
****************************************************/
   y=Mat<Type> (out_rows,out_cols);
   Mat<Type> yi(nrows,out_cols);
   
   Mat<Type> tempr(1,ncols*Urow+1);
   
   for(i=0 ; i < nrows ; i++){
      tempr.fill((Type)0);
/*************************
      tempr = Zeros<Type>(1,ncols*Urow + 1);
**************************/
      tempr(0) = _el[i][0];
      for(n=1 ; n < ncols ; n++){
	 tempr(n*Urow) = _el[i][n];
	 for(m=1 ; m < Urow ; m++){
	    tempr((n-1)*Urow + m) = 
	       tempr((n-1)*Urow) + (tempr(n*Urow)-tempr((n-1)*Urow))*m/Urow;
	 }
      }
      for(n=0 ; n < out_cols ; n++)
	 yi(i,n) = tempr(n*Drow);
   }
   
   Mat<Type> tempc(1,nrows*Ucol + 1);
   
   for(i=0 ; i < out_cols ; i++){
/***************changed 6/23/94****************   
      tempc = Zeros<Type>(1,nrows*Ucol + 1);
**********************************************/
      tempc.fill((Type)0);
      tempc(0) = yi(0,i);
      for(n=1 ; n < nrows ; n++){
	 tempc(n*Ucol) = yi(n, i);
	 for(m=1 ; m < Ucol ; m++){
	    tempc((n-1)*Ucol + m) = 
	       tempc((n-1)*Ucol) + (tempc(n*Ucol)-tempc((n-1)*Ucol))*m/Ucol;
	 }
      }
      for(n=0 ; n < out_rows ; n++)
	 y(n,i) = tempc(n*Dcol);
   }
}

//
//-------------------------//
//

template <class Type>
Mat<Type> 
Mat<Type>::filter(Mat<Type>& B,Mat<Type>& A) const
{
   unsigned     n, i, P, Q, N;
   Type  temp;
   
   if ((B._rows != 1) && (B._cols != 1)){
      cerr << "error in filter: numerator not a vector." << endl;
      cerr << B._rows << "x" << B._cols << endl;
      exit(1);
   }
   if ((A._rows != 1) && (A._cols != 1)){
      cerr << "error in filter: denominator not a vector." << endl;
      cerr << A._rows << "x" << A._cols << endl;
      exit(1);
   }
   if ((_rows != 1) && (_cols != 1)){
      cerr << "error in filter: input not a vector." << endl;
      cerr << _rows << "x" << _cols << endl;
      exit(1);
   }
   if (A(0) == 0){
      cerr << "error in filter: first element in the denominator" << endl;
      cerr << "      must be non-zero." << endl;
      exit(1);
   }
   P = A.length();
   P=P-1;   
   Q = B.length();   
   Q=Q-1;
   N = this->length();
   temp = A(0);
   if (A(0) != 1.0){
      for (i=0 ; i <= Q ; i++)
	 B(i) = B(i)/temp;
      for (i=0 ; i <= P ; i++)
	 A(i) = A(i)/temp;
   }
   
   Mat<Type> y(_rows,_cols);
   Mat<Type> pastx(1,Q+1);   
   Mat<Type> pasty(1,P+1);
   
   for(n=0 ; n < N ; n++){
      y(n) = 0;
      pastx(0) = this->operator () (n);
      for (i=1 ; i <= P ; i++)
	 y(n) = y(n) - A(i)*pasty(i);
      for (i=0 ; i <= Q ; i++)
	 y(n) = y(n) + B(i)*pastx(i);
      pasty(0) = y(n);
      for(i=P ; i > 0 ; i--)
	 pasty(i) = pasty(i-1);
      for(i=Q ; i > 0 ; i--)
	 pastx(i) = pastx(i-1);
   }
   
   return(y);
}

/***************************morphology functions*******************************/
#ifdef USE_DBLMAT
template <class Type>
Mat<Type>
Mat<Type>::erode(const Mat<double>& strel) const
{
  unsigned strelHeight = strel.getrows();
  unsigned strelWidth  = strel.getcols();

  if (((strelHeight == 1) && (strelWidth == 1)) || !strelHeight || !strelWidth)
    return Mat<Type>(*this);

  Mat<Type> padMatrix(pad(strelHeight/2, strelWidth/2));
  Mat<Type> result(_rows, _cols);

  Type   *padMatrixPtr1 = padMatrix._el[0];
  const double *strelStart = strel.getEl()[0];
  Type   *resultPtr     = result._el[0];

  unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
  unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
     
  for (unsigned y = _rows; y != 0; y--) {
    for (unsigned x = _cols; x != 0; x--) {
      Type   *padMatrixPtr2 = padMatrixPtr1;
      double  minimum = MAXDOUBLE;
      const double *strelPtr = strelStart;
      for (register unsigned j = strelHeight; j != 0; j--) {
	for (register unsigned i = strelWidth; i != 0; i--) {
	  if (*strelPtr >= 0)
	    minimum = MIN(minimum, (double) *padMatrixPtr2 - *strelPtr);
	  strelPtr++;
	  padMatrixPtr2++;
	}
	padMatrixPtr2 += padMatrixPtr2incr;
      }
      *resultPtr = Type(minimum);
      
      resultPtr++;
      padMatrixPtr1++;
    }
    padMatrixPtr1 += padMatrixPtr1incr;
  }
  
  return result;
}

template <class Type>
Mat<Type>
Mat<Type>::dilate(const Mat<double>& strel) const
{
   unsigned strelHeight = strel.getrows();
   unsigned strelWidth  = strel.getcols();

   if (((strelHeight == 1) && (strelWidth == 1)) || !strelHeight || !strelWidth)
      return Mat<Type>(*this);

   unsigned r=0;
   unsigned c=0;

   if ((strelHeight % 2) == 0) {
      strelHeight++;
      r++;
   }
   
   if ((strelWidth % 2) == 0) {
      strelWidth++;
      c++;
   }
   
   Mat<double> newStrel(strelHeight, strelWidth, -1);
   newStrel.insert(strel.rotate180(), r, c);
   
   Mat<Type> padMatrix(pad(strelHeight/2, strelWidth/2));
   Mat<Type> result(_rows, _cols);
   
   Type   *padMatrixPtr1 = padMatrix._el[0];
   const double *strelStart    = newStrel.getEl()[0];
   Type   *resultPtr     = result._el[0];
   
   unsigned padMatrixPtr2incr = padMatrix._cols - strelWidth;
   unsigned padMatrixPtr1incr = (strelWidth/2) * 2;
   
   for (unsigned y = _rows; y != 0; y--) {
      for (unsigned x = _cols; x != 0; x--) {
	 Type   *padMatrixPtr2 = padMatrixPtr1;
	 double  maximum = -MAXDOUBLE;
	 const double *strelPtr = strelStart;
	 for (register unsigned j = strelHeight; j != 0; j--) {
	    for (register unsigned i = strelWidth; i != 0; i--) {
	       if (*strelPtr >= 0)
		  maximum = MAX(maximum, (double) *padMatrixPtr2 + *strelPtr);
	       strelPtr++;
	       padMatrixPtr2++;
	    }
	    padMatrixPtr2 += padMatrixPtr2incr;
	 }
	 *resultPtr = Type(maximum);
	 
	 resultPtr++;
	 padMatrixPtr1++;
      }
      padMatrixPtr1 += padMatrixPtr1incr;
   }
   
   return result;
}

#endif // USE_DBLMAT

/******************************Old code***************************************
template <class Type>
Mat<Type> 
Mat<Type>::filter(Mat<Type>& B, Mat<Type>& A, Mat<Type>& x)
{
   unsigned     n, i, P, Q, N;
   Type  temp;
   
   if ((B._rows != 1) && (B._cols != 1)){
      cerr << "error in filter: numerator not a vector." << endl;
      cerr << B._rows << "x" << B._cols << endl;
      exit(1);
   }
   if ((A._rows != 1) && (A._cols != 1)){
      cerr << "error in filter: denominator not a vector." << endl;
      cerr << A._rows << "x" << A._cols << endl;
      exit(1);
   }
   if ((x._rows != 1) && (x._cols != 1)){
      cerr << "error in filter: input not a vector." << endl;
      cerr << x._rows << "x" << x._cols << endl;
      exit(1);
   }
   if (A(0) == 0){
      cerr << "error in filter: first element in the denominator" << endl;
      cerr << "      must be non-zero." << endl;
      exit(1);
   }
   P = length(A) - 1;   Q = length(B) - 1;   N = length(x);
   temp = A(0);
   if (A(0) != 1.0){
      for (i=0 ; i <= Q ; i++)
	 B(i) = B(i)/temp;
      for (i=0 ; i <= P ; i++)
	 A(i) = A(i)/temp;
   }
   
   Mat<Type> y(x._rows, x._cols);
   Mat<Type> pastx(1,Q+1);   
   Mat<Type> pasty(1,P+1);   
   pastx = Zeros<Type>(1,Q+1);
   pasty = Zeros<Type>(1,P+1);
   
   for(n=0 ; n < N ; n++){
      y(n) = 0;
      pastx(0) = x(n);
      for (i=1 ; i <= P ; i++)
	 y(n) = y(n) - A(i)*pasty(i);
      for (i=0 ; i <= Q ; i++)
	 y(n) = y(n) + B(i)*pastx(i);
      pasty(0) = y(n);
      for(i=P ; i > 0 ; i--)
	 pasty(i) = pasty(i-1);
      for(i=Q ; i > 0 ; i--)
	 pastx(i) = pastx(i-1);
   }
   
   return(y);
}
****************************************************************************/



//Need to be modified 
//
//-------------------------//
//

/****************************************************************
template <class Type>
void 
Mat<Type>::radialscale(Mat<Type>& scale, Mat<Type>& y) const
{  //To recover old code look in rmold.c
   unsigned    slen, out_rows, out_cols, nrows, ncols;
   unsigned    i, j;
   Type smax, smin, fPI;
   Type xrc, xcc, yrc, ycc;
   Type xang, yang, xrad, yrad, xr, xc;
   unsigned rloc,cloc;
   
   fPI = 4.0*atan(1.0);
   slen = 360;
   nrows = _rows;   ncols = _cols;
   smin = min(scale);     smax = max(scale);
   out_rows = (unsigned)(floor(smax*nrows + 0.5));
   out_cols = (unsigned)(floor(smax*ncols + 0.5));
   y=Mat<Type>(out_rows,out_cols);
   
   xrc = (Type)(nrows)/2.0;     xcc = (Type)(ncols)/2.0;
   yrc = (Type)(out_rows)/2.0;   ycc = (Type)(out_cols)/2.0;
   
   for(i=0 ; i < out_rows ; i++)
      for(j=0 ; j < out_cols ; j++){
	 yrad = (Type)(::sqrt((i-yrc)*(i-yrc) + (j-ycc)*(j-ycc)));
	 yang = (Type)floor(atan2((j-ycc) , (i-yrc))*180/fPI);
	 if (yang < (Type)0) 
	    yang = yang + (Type)360;
	 xrad = yrad/scale((int)(yang));
	 xang = yang;
	 xr   = (Type)floor(xrc + xrad*(Type)cos(xang*fPI/(Type)180) );
	 xc   = (Type)floor(xcc + xrad*(Type)sin(xang*fPI/(Type)180));
	 if ((xr<(Type)0) || (xr>=(Type)nrows) || (xc<(Type)0) || (xc>=(Type)ncols))
	    y(i,j) = (Type)0;
	 else
	    y(i,j) = _el[(unsigned)xr][(unsigned)xc];
	 
      }
}
******************************************************************************/

//
//-------------------------//
//
//
//------------------------/
//

#ifdef HAVE_MATLAB
/*****************************************************************************/
//This function goes to the named mat file and searches for the given variable 
//Next, if possible, it reinitializes the calling object to represent the
//same matrix as that defined by the Matlab variable.  The Matlab matrix
//is stored in double format.  The function Matlab2Mat will convert the 
//Matlab data to the type of its calling object (this).  This may mean a loss
//of precision.
template <class Type>
Boolean
Mat<Type>::loadMatlab(const char *filename, const char *varname)
{
  if (!filename || !varname) {
    cerr << "Mat<Type>::loadMatlab(): no filename or variable name specified" << endl;
    return FALSE;
  }

  //Check to ensure that the file can be opened successfully for reading
  MATFile *infile = matOpen((char *) filename,"r");
  if(infile==NULL) {
    cerr << "Mat<Type>::loadMatlab(): Can't open file: " << filename << endl;
    return FALSE;
  }
   
  //Assign pointer to the Matrix object allocated from the file and ensure that
  //the named variable was found
  Matrix *MatlabMatrix = matGetMatrix(infile, (char *) varname);
  if(MatlabMatrix==NULL) {
    cerr << "Mat<Type>::loadMatlab(): Could not find: " << varname << " in file: " 
	 << filename << endl;
    //close the file and quit 
    matClose(infile);
    return FALSE;
  }
  
  //Check to see if the Matlab Matrix was full or sparse, dcomplex or real
  //if (mxIsDcomplex(MatlabMatrix) == 0)
  if (mxIsComplex(MatlabMatrix) == 0)
    *this = matlab2Mat(MatlabMatrix);
  else {
    cerr<<"Mat<Type>::loadMatlab(): The variable: "<<varname<<" in the file: "<<filename;
    cerr<<" is dcomplex.  Create a dcomplex Mat object and call its";
    cerr<<" loadMatlab function."<<endl;
    matClose(infile);
    mxFreeMatrix(MatlabMatrix);
    return FALSE;
  }
  
  //We were successful
  matClose(infile);
  mxFreeMatrix(MatlabMatrix);

  return TRUE;
}

#endif // HAVE_MATLAB
//
//-------------------------//
//
template <class Type>
Boolean 
Mat<Type>::loadRaw(const char *filename, unsigned nrows, unsigned ncols)
{
  InputFile matrixFile(filename);

  if (!matrixFile) {
    cerr << "Error in loadRaw: error opening file." << endl;
    return FALSE;
  }

  _checkMatrixDimensions(filename, nrows, ncols);
  
  if ((nrows && (nrows != _rows)) || (_cols && (ncols != _cols))) {
    _rows = _maxrows = nrows;
    _cols = _maxcols = ncols;
    _allocateEl();
  }
  
  if (!matrixFile.stream().read((char *)_el + _maxrows*sizeof(Type *),
				_maxrows*_maxcols*sizeof(Type)))
    return FALSE;
  
  return TRUE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::loadAscii(const char *filename)
{
  InputFile infile(filename);
  if (!infile){
    cerr << "Error in loadAsccii: error opening file." << endl;
    return FALSE;
  }
  
  if (!(infile >> _rows >> _cols))
    return FALSE;
  
  _maxrows = _rows;
  _maxcols = _cols;
  
  _allocateEl();
  
  for(unsigned i=0; i < _rows; i++) {
    for(unsigned j=0; j < _cols ; j++)
      if (!(infile >> _el[i][j]))
	return FALSE;
  }
  
  return TRUE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::load(const char *filename, int type, const char *varname)
{
  Boolean status = FALSE;

  switch(type) {
  case RAW:
    status = loadRaw(filename);
    break;
  case ASCII:
    status = loadAscii(filename);
    break;
  case MATLAB:
#ifdef HAVE_MATLAB
    status = loadMatlab(filename, varname);
#else
    cerr << "Can't load MATLAB files. Please compile with -DHAVE_MATLAB" << endl;
#endif
    break;
  default:
    cerr << "Unrecognized type for loading" << endl;
    break;
  }

  return status;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::save(const char *filename, int type, const char *varname) const
{
  Boolean result;		// bert - added

  switch(type) {
  case RAW:
    result = saveRaw(filename);
    break;
  case ASCII:
    result = saveAscii(filename);
    break;
  case MATLAB:
#ifdef HAVE_MATLAB
    result = saveMatlab(filename, varname);
#else
    cerr << "Can't save MATLAB files. Please compile with -DHAVE_MATLAB" << endl;
#endif
    break;
  default:
    cerr << "Unrecognized type for saving" << endl;
    result = 0;
    break;
  }
  return (result);		// bert - return a value.
}

//
//-------------------------//
//
#ifdef HAVE_MATLAB
/*****************************************************************************/
// Saves the calling object as a Matlab readable matrix to the named file with the
// given variable name. Since Matlab stores in double format, the calling object is 
// first converted to a temporary double array. The third argument lets the user
// update a file("u") or rewrite an existing file or create a new one("w")
template <class Type>
Boolean
Mat<Type>::saveMatlab(const char *filename, const char *varname, const char *option) const
{
  if (!_rows || !_cols) {
    cerr << "Attempt to save empty matrix as matlab file. Not saved" << endl;
    return FALSE;
  }

  //Create a temporary double array
  double *temp = 0;
  allocateArray(_rows*_cols, temp);
  if (!temp) {
    cerr << "Couldn't allocate temp array; matrix not saved" << endl;
    return FALSE;
  }

  // Fill temp array columnwise (this is Matlab's storage convention)
  double *tempPtr = temp;
  for (unsigned col = 0; col < _cols; col++)
    for (unsigned row = 0; row < _rows; row++)
      *tempPtr++ = _el[row][col];
  
  // Write the matrix out
  Boolean status = ::saveMatlab(filename, varname, option, _rows, _cols, temp);

  // Free temporary array
  freeArray(temp);

  return status;
}
#endif // HAVE_MATLAB

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::saveRaw(const char *filename) const
{
   ofstream outfile(filename);
   if (!outfile){
      cerr << "Error in saveRaw: error opening file." << endl;
      return FALSE;
   }
   
//outfile.write((unsigned char *) &nrows, sizeof(unsigned));
//outfile.write((unsigned char *) &ncols, sizeof(unsigned));
   
   outfile.write((char *)_el + _maxrows*sizeof(Type *) , _maxrows*_maxcols*sizeof(Type));
   
   outfile.close();

   return (outfile) ? TRUE : FALSE;
}

//
//-------------------------//
//
template <class Type>
Boolean
Mat<Type>::saveAscii(const char *filename) const
{
   unsigned   i, j;
   
   ofstream outfile(filename);
   if (!outfile){
      cerr << "Error in saveAscciifile: error opening file." << endl;
      return FALSE;
   }
   
   outfile << _rows << " " << _cols << endl;
   
   for(i=0; i < _rows; i++)
   {
      for(j=0; j < _cols ; j++)
	 outfile << _el[i][j] << " ";
      outfile << endl;
   }
   
   outfile.close();

   return (outfile) ? TRUE : FALSE;
}

//
//-------------------------//
//
 

//
//-------------------------//
//

//
// Private functions start
//
template <class Type>
void
Mat<Type>::_allocateEl()
{
  if (_el) {
#ifdef DEBUG
    cout << "Freeing allocated memory at " << _el << endl;
#endif
    if( _el[0] ) {
      delete [] _el[0];  // delete data block
    }
    delete [] _el; // delete row of pointers  
  }
  _el = 0;

  unsigned int nBytes = _maxrows*_maxcols*sizeof(Type);

  if (nBytes) {
    typedef Type * TypePtr;
    _el = (Type **) new TypePtr[_maxrows];
    assert(_el);
    _el[0] = (Type *) new Type[_maxrows*_maxcols];
    assert(_el[0]);
    memset(_el[0], 0, nBytes);  //set all elements to zero

#ifdef DEBUG
    cout << "Allocated " << nBytes << " bytes at " <<  _el << endl;
#endif
   
    for (unsigned i=1 ; i < _maxrows ; i++)
      _el[i] = (Type *)(&(_el[i-1][_maxcols]));
  }
}

/*************************************************************************
This function will actually copy another matrix to a subsection of the
calling matrix
Since it modifies data internally, it has been declared private
**************************************************************************/

template <class Type> 
void 
Mat<Type>::_modify_sub_section(unsigned r1,unsigned r2,unsigned c1,unsigned c2,const Mat<Type>& B)
{
  if ((r1 > r2) || (c1 > c2) || (r2 >= _rows) || (c2>= _cols)) { // (bert)
      cerr << "Error in cropting: improper row or column sizes." << endl;
      cerr << r1 << " to " << r2 << " and" << endl;
      cerr << c1 << " to " << c2 << endl;
      exit(1);
   }
//If the selected subsection and B object are not the same size, then done   
   if((B.getrows() != (r2-r1+1)) || (B.getcols() != (c2-c1+1)))
   {
     cerr<<"Error:Input Matrix and subsection selections don't argree in size";
     cerr<<endl;
   } 
   
   
   for (unsigned i=r1; i<=r2 ; i++)
     for (unsigned j=c1; j<=c2; j++)
	   _el[i][j] = B(i-r1,j-c1);
          
}

template <class Type>
void
Mat<Type>::_checkMatrixDimensions(const char *path, unsigned& nrows,
				  unsigned& ncols) const
{
  if (Path(path).hasCompressedExtension())
    return;

  struct stat buf;
  stat(path, &buf);
  inferDimensions(buf.st_size/sizeof(Type), nrows, ncols);
}

template <class Type>
Mat<Type>&
Mat<Type>::_fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc)
{
  using std::max;		// (bert) force it to use the right max()
  using std::abs;		// (bert) and force the right abs()!

  // Verify dimensions of FFT
  if ((nrows > 1) && ((nrows < _rows) || !isPowerOf2(nrows))) {
    cerr << "Warning! Mat<Type>::fft():" << endl
	 << "  Requested # _rows for FFT (" << nrows << ") invalid;" << endl;
    nrows = max(unsigned(4), nextPowerOf2(max(nrows, _rows)));
    cerr << "  increased to " << nrows << endl;
  }

  if ((ncols > 1) && ((ncols < _cols) || !isPowerOf2(ncols))) {
    cerr << "Warning! Mat<Type>::fft():" << endl
	 << "  Requested # _cols for FFT (" << ncols << ") invalid;" << endl;
    ncols = max(unsigned(4), nextPowerOf2(max(ncols, _cols)));
    cerr << "  increased to " << ncols << endl;
  }

  Boolean doX = (ncols != 1) && (_cols > 1) ? TRUE : FALSE;
  Boolean doY = (nrows != 1) && (_rows > 1) ? TRUE : FALSE;

  if (!nrows)
    nrows = max(unsigned(4), nextPowerOf2(_rows));
  else if (nrows == 1)
    nrows = _rows;

  if (!ncols)
    ncols = max(unsigned(4), nextPowerOf2(_cols));
  else if (ncols == 1)
    ncols = _cols;

  // Pad matrix to the final FFT dimensions
  pad(nrows, ncols, (nrows - _rows)/2, (ncols - _cols)/2, 0);

  // Allocate temporary arrays
  double *real = 0;
  double *imag = 0;
  unsigned maxDim = max(doX ? _cols : 1, doY ? _rows : 1);

  allocateArray(maxDim, real);
  allocateArray(maxDim, imag);
  assert(real && imag);

  unsigned row, col;

  // Take 1D FFT in X (row) direction
  if (doX)
    for(row = 0; row < _rows; row++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      Type *sourcePtr = _el[row];
      for(col = _cols; col != 0; col--) {
	dcomplex value(*sourcePtr++);
	*realPtr++ = value.real();
	*imagPtr++ = value.imag();
      }

      // Calculate 1D FFT
      fftFunc(_cols, real, imag);

      // Put results back
      realPtr   = real;
      imagPtr   = imag;
      sourcePtr = _el[row];
      for(col = _cols; col != 0; col--)
	*sourcePtr++ = Type(abs(dcomplex(*realPtr++, *imagPtr++)));
    }

  // Take 1D FFT in Y (column) direction
  if (doY)
    for(col = 0; col < _cols; col++) {
      // Fill temporary FFT array
      double  *realPtr   = real;
      double  *imagPtr   = imag;
      Type    *sourcePtr = _el[0] + col;
      for(row = _rows; row != 0; row--) {
	dcomplex value(*sourcePtr);
	*realPtr++ = value.real();
	*imagPtr++ = value.imag();
	sourcePtr += _cols;
      }

      // Calculate 1D FFT
      fftFunc(_rows, real, imag);

      // Put results back
      realPtr   = real;
      imagPtr   = imag;
      sourcePtr = _el[0] + col;
      for(row = _rows; row != 0; row--) {
	*sourcePtr = Type(abs(dcomplex(*realPtr++, *imagPtr++)));
	sourcePtr += _cols;
      }
    }

  // Free temporary arrays
  freeArray(real);
  freeArray(imag);

  return *this;
}

template <class Type>
Mat<Type>&
Mat<Type>::ifft(unsigned nrows, unsigned ncols)
{
  _fft(nrows, ncols, ::ifft);
  unsigned factor = 1;
  if (nrows != 1) factor *= _rows;
  if (ncols != 1) factor *= _cols;
  
  return *this /= factor;
}

#ifdef HAVE_MATLAB
//This function converts a Matlab matrix object to a Mat object
//The first argument is a pointer to a Matlab matrix object
//This function checks the storage format of the Matlab object
//(Full or Sparse) and does the appropriate conversions to return a Mat
//object of the calling object's type.
template <class Type>
Mat<Type>
Mat<Type>::matlab2Mat(Matrix *MatlabMatrix) const
{  
   //Make sure that pointer is not NULL
   if(MatlabMatrix==NULL)
     {
       cerr<<"The input argument points to a non existent object(NULL)"<<endl;
       exit(1);
     }
   
   //See if the object is real or dcomplex
   //if(mxIsDcomplex(MatlabMatrix))
   if(mxIsComplex(MatlabMatrix))
     {
       cerr<<"The Matlab matrix object is dcomplex.  Create a dcomplex Mat ";
       cerr<<"object first then repeat the call"<<endl;
       exit(1);
     }
   
   unsigned Matlab_rows=mxGetM(MatlabMatrix);
   unsigned Matlab_cols=mxGetN(MatlabMatrix);
   Mat<Type> Temp(Matlab_rows,Matlab_cols);
   double *MatlabDATA=mxGetPr(MatlabMatrix);
   //See if the object is Full or Sparse and do corresponding action
   if(mxIsFull(MatlabMatrix))
     { //Since the Matlab storage representation of the object
       //is exactly the transpose of the Mat representation, The object
       //is initialized with its indices reversed.
       
       for(unsigned i=0;i<Matlab_cols;i++)
          for(unsigned j=0;j<Matlab_rows;j++)
            {
             Temp._el[j][i]= (Type) *MatlabDATA;
             MatlabDATA++;
            } 
       return(Temp);
        
     }
   else if(mxIsSparse(MatlabMatrix))
     {
       //This code was derived from code form the external interface guide for 
       // Matlab pg. 7 (Copyright Math Works)
       int *ir=mxGetIr(MatlabMatrix);
       int *jc=mxGetJc(MatlabMatrix);
       int temp_index;
       for(int j=0; j<Matlab_cols;j++)
        {
          temp_index=jc[j+1];
          for(int k=jc[j];k<temp_index;k++)
            Temp._el[ir[k]][j] = (Type) MatlabDATA[k];
        }
       
       return(Temp);
     }

   return(NULL);
}

template <class Type>
Matrix *
Mat<Type>::mat2Matlab(const char *name) const
{ 
  Matrix *CopyofThis;

  //Allocate the memory for the Matlab matrix and set its name
  CopyofThis=mxCreateFull(_rows,_cols,0);
  mxSetName(CopyofThis,name);

  double *MatlabData=mxGetPr(CopyofThis);
  for(unsigned i=0;i< _cols;i++)
    for(unsigned j=0;j<_rows;j++)
     {
      *MatlabData= (double) _el[j][i];
       MatlabData++;
     } 
  
  return(CopyofThis);
}

#endif // HAVE_MATLAB

//
// Private functions end
//

template <class Type>
ostream&
operator << (ostream& os, const Mat<Type>& A)
{
  return A.display(os);
}

// 
// Type conversions
//

#ifdef USE_DBLMAT
template <class Type> Mat<double> asDblMat(const Mat<Type>& A)
{
  Mat<double> cast(A.getrows(), A.getcols());
  double     *castPtr = (double *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++ = asDouble(*aPtr++);
  
  return cast;
}
#endif // USE_DBLMAT

//
//-------------------------// 
//

#ifdef USE_FLMAT
template <class Type> Mat<float> asFlMat(const Mat<Type>& A)
{
  Mat<float>  cast(A.getrows(), A.getcols());
  float      *castPtr = (float *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++ = (float) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_INTMAT
template <class Type> Mat<int> asIntMat(const Mat<Type>& A)
{
  Mat<int>    cast(A.getrows(), A.getcols());
  int        *castPtr = (int *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (int) asDouble(*aPtr++);

  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_UINTMAT
template <class Type> Mat<unsigned int> asUIntMat(const Mat<Type>& A)
{
  Mat<unsigned int> cast(A.getrows(), A.getcols());
  unsigned        *castPtr = (unsigned *) cast.getEl()[0];
  const Type      *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_SHMAT
template <class Type> Mat<short> asShMat(const Mat<Type>& A)
{
  Mat<short>  cast(A.getrows(), A.getcols());
  short      *castPtr = (short *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (short) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_USHMAT
template <class Type> Mat<unsigned short> asUShMat(const Mat<Type>& A)
{
  Mat<unsigned short>  cast(A.getrows(), A.getcols());
  unsigned short      *castPtr = (unsigned short *) cast.getEl()[0];
  const Type          *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned short) asDouble(*aPtr++);
  
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_CHRMAT
template <class Type> Mat<char> asChrMat(const Mat<Type>& A)
{
  Mat<char>   cast(A.getrows(), A.getcols());
  char       *castPtr = (char *) cast.getEl()[0];
  const Type *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (char) asDouble(*aPtr++);
   
  return cast;
}
#endif

//
//-------------------------// 
//

#ifdef USE_UCHRMAT
template <class Type> Mat<unsigned char> asUChrMat(const Mat<Type>& A)
{
  Mat<unsigned char> cast(A.getrows(), A.getcols());
  unsigned char     *castPtr = (unsigned char *) cast.getEl()[0];
  const Type        *aPtr    = A.getEl()[0];
  for (unsigned i=A.nElements(); i; i--)
    *castPtr++  = (unsigned char) asDouble(*aPtr++);
  
  return cast;
}
#endif

// (bert) added this declaration here to make the compiler happier.  This
// way it doesn't see the default parameters as being listed as part of
// the function definition.  I think this was causing a compilation problem
// Jason found.
//
template<class Type>
SimpleArray<Type> array(const Mat<Type>& A, Type minVal, Type maxVal)
{
    return A.array(minVal, maxVal); 
}

//////////////////////////////////////////////////////////////////////

#endif
