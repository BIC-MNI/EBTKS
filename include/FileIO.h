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
$RCSfile: FileIO.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef FILE_IO_H
#define FILE_IO_H

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include "Path.h"
//#include "popen.h"

#ifdef HAVE_MATLAB
extern "C"{    
  #include"mat.h" 
}
#endif

class InputFile {
  //ipopen *_ipipe;
  istream *_ipipe;
  
public:
// Constructors/destructor
  InputFile();
  InputFile(const Path& path) { _ipipe = 0; attach(path); }
  InputFile(istream *pipe)   { _ipipe = pipe; }
  ~InputFile()                { close(); }

  operator void *() const;
  operator istream& () { return *_ipipe; }

// Set functions
  Boolean  attach(const Path&);
  istream& skip(unsigned nBytes);
  Boolean  close();

// Get functions
  istream& stream()     { return *_ipipe; }
  Boolean operator ! () { return (!_ipipe || !*_ipipe); }

// Input operators
  istream& operator >> (char *val) { return *_ipipe >> val; }
  istream& operator >> (char& val) { return *_ipipe >> val; }
  istream& operator >> (short& val) { return *_ipipe >> val; }
  istream& operator >> (int& val) { return *_ipipe >> val; }
  istream& operator >> (long& val) { return *_ipipe >> val; }
  istream& operator >> (float& val) { return *_ipipe >> val; }
  istream& operator >> (double& val) { return *_ipipe >> val; }
  istream& operator >> (unsigned char *val) { return *_ipipe >> val; }
  istream& operator >> (unsigned char& val) { return *_ipipe >> val; }
  istream& operator >> (unsigned short& val) { return *_ipipe >> val; }
  istream& operator >> (unsigned int& val) { return *_ipipe >> val; }
  istream& operator >> (unsigned long& val) { return *_ipipe >> val; }
  istream& operator >> (streambuf *val) { return *_ipipe >> val; }
  istream& operator >> (istream& (*func)(istream&)) { return *_ipipe >> func; }
  istream& operator >> (ios& (*func)(ios&)) { return *_ipipe >> func; }
};

/******************
 * OutputFile class
 ******************/

class OutputFile : public ofstream {
  Path _path;
  int  _compress;

public:
  static const Boolean NO_COMPRESS, COMPRESS;

  OutputFile(const Path&, int mode = ios::out, int compress = COMPRESS);
  ~OutputFile();

// Set functions
  void compress(Boolean);
};

/***************************
 * C file handling functions
 ***************************/

FILE    *openFile(const Path&, const char *);
void     closeFile();

#ifdef HAVE_MATLAB
MATFile *openMatlabFile(const Path&, const char *);
void     closeMatlabFile(MATFile *);
#endif

#endif
