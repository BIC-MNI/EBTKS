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
$RCSfile: test.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:24 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "CachedArray.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

void printContents(const SimpleArray<float>& A);

/*
 * Little utility functions not necessary for median-finding.
 */

float current_cpu_usage (void)
{
   struct rusage  usage;

   getrusage (RUSAGE_SELF, &usage);
   return (float) usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec / 1e6);
}

float prev_tic;

void tic (char *msg)
{
   if (msg != NULL)
   {
     cout << "timing operation: " << msg << endl;
   }
      
   prev_tic = current_cpu_usage ();
}

void toc (void)
{
  cout << "time elapsed = " << current_cpu_usage() - prev_tic << " sec" << endl;
}

void
main(int argc, char **argv) {
  unsigned size = 20;
  unsigned nBlocks = 2;
  unsigned blockSize = 7;

  if (argc > 3)
    blockSize = atoi(argv[3]);
  if (argc > 2)
    nBlocks = atoi(argv[2]);
  if (argc > 1)
    size = atoi(argv[1]);

  SimpleArray<float> *array = new CachedArray<float>(size, nBlocks, blockSize);
  //SimpleArray<float> *array = new SimpleArray<float>(size);
  assert(array);

  SimpleArray<float>& A = *array;

  for (unsigned i = 0; i < size; i++)
    A[i] = SQR(gauss(10, 10));

  printContents(A);

  CachedArray<float> B;

  fstream stream("test.ar", ios::in|ios::out);
  A.saveBinary(stream);
  stream.seekg(0);
  B.loadBinary(stream, size);

  cout << "A: " << A << endl;
  cout << "B: " << B << endl;

  stream.seekg(0);
  A.saveAscii(stream);
  stream.seekg(0);
  B.clear(9);
  B.loadAscii(stream, size);

  cout << "A: " << A << endl;
  cout << "B: " << B << endl;

  A.clear(10);

  printContents(A);

  A.newSize(10);

  A.newSize(50);
  A.clear(8);

  B = A;

  A.clear(15);

  B += A;

  printContents(B);

  B /= 0.4;
  B += 3;
  B -= 3;
  B *= 0.4;

  printContents(B);
}

void printContents(const SimpleArray<float>& A)
{
  //  cout << endl << "Contents: " << A << endl;

  cout << endl << "Size: " << A.size() << endl;
  unsigned i;
  cout << "Min: " << A.min(&i) << " at " << i << endl;
  cout << "Max: " << A.max(&i) << " at " << i << endl;
  cout << "Sum: " << A.sum() << endl;
  cout << "Sum2: " << A.sumSqr() << endl;
  cout << "Mean: " << A.mean() << endl;
  cout << "Prod: " << A.prod() << endl;
  cout << "Prod2: " << A.prodSqr() << endl;
  cout << "Var: " << A.var() << endl;
  //A.resetHitRate();
  cout << "Std: " << A.std() << endl; //" HR: " << A.hitRate() << endl;
  //A.resetHitRate();
  tic("median");
  cout << "Median: " << A.median() << endl; //" HR: " << A.hitRate() << endl;
  toc();
}




