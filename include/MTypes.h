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
$RCSfile: MTypes.h,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2003-04-16 18:00:29 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef MTYPES_H
#define MTYPES_H

#include <values.h>

typedef unsigned char  Byte;
//typedef unsigned short Grey;
//typedef unsigned long  Full;
typedef short Grey;
typedef long  Full;

typedef char  Boolean; // In accordance with X's definition

const int BYTE_MIN = 0;
const int BYTE_MAX = 255;
const int GREY_MIN = -MAXSHORT;
const int GREY_MAX = MAXSHORT;

inline Byte 
clipByte(double value) { 
  Byte newValue = (Byte) value; 
  if (value < (double) BYTE_MIN) 
    newValue = BYTE_MIN;
  else if (value > (double) BYTE_MAX) 
    newValue = BYTE_MAX;
  return newValue;
}

inline Byte 
clipByte(Full value) { 
  Byte newValue = (Byte) value; 
  if (value < (double) BYTE_MIN) 
    newValue = BYTE_MIN;
  else if (value > (double) BYTE_MAX) 
    newValue = BYTE_MAX;
  return newValue;
}

inline Byte 
clipByte(int value) { 
  Byte newValue = (Byte) value; 
  if (value < (double) BYTE_MIN) 
    newValue = BYTE_MIN;
  else if (value > (double) BYTE_MAX) 
    newValue = BYTE_MAX;
  return newValue;
}

inline Grey 
clipGrey(double value) { 
  Grey newValue = (Grey) value; 
  if (value < (double) GREY_MIN) 
    newValue = GREY_MIN;
  else if (value > (double) GREY_MAX) 
    newValue = GREY_MAX;
  return newValue;
}

inline Grey 
clipGrey(Full value) { 
  Grey newValue = (Grey) value; 
  if (value < (double) GREY_MIN) 
    newValue = GREY_MIN;
  else if (value > (double) GREY_MAX) 
    newValue = GREY_MAX;
  return newValue;
}

inline Grey 
clipGrey(int value) { 
  Grey newValue = (Grey) value; 
  if (value < (double) GREY_MIN) 
    newValue = GREY_MIN;
  else if (value > (double) GREY_MAX) 
    newValue = GREY_MAX;
  return newValue;
}

#endif
