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
$RCSfile: OpTimer.cc,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "OpTimer.h"
#include "trivials.h"
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>

const char *OpTimer::_TIME_STRINGS[] = { "CPU", "SYS", "USR" };

const int   OpTimer::CPU = 0;
const int   OpTimer::SYS = 1;
const int   OpTimer::USR = 2;

OpTimer::OpTimer(int timeType, const char *operation, unsigned N, 
		 unsigned reportInterval) 
{
  _os = &cout;
  this->timeType(timeType);
  _operation = 0;
  _verbose = TRUE;
  _newOperation(operation);
  tic(N, reportInterval);
}

OpTimer::OpTimer(const char *operation, unsigned N, unsigned reportInterval) 
{
  _os = &cout;
  timeType(USR);
  _operation = 0;
  _verbose = TRUE;
  _newOperation(operation);
  tic(N, reportInterval);
}

OpTimer::~OpTimer()
{
  if (!_N)
    toc();

  if (_verbose)
    *_os << flush;
}

double 
OpTimer::operator () (int timeType, const char *operation, unsigned N, 
		      unsigned reportInterval)
{
  this->timeType(timeType);
  _newOperation(operation);
  return tic(N, reportInterval);
}

double 
OpTimer::operator () (const char *operation, unsigned N, unsigned reportInterval)
{
  _newOperation(operation);
  return tic(N, reportInterval);
}

void
OpTimer::timeType(int time)
{
  _timeType = time;

  switch (time) {
  case CPU:
    timeFunction(_CPUtime);
    break;
  case SYS:
    timeFunction(_SYStime);
    break;
  case USR:
    timeFunction(_USRtime);
    break;
  default:
    cerr << "Warning! Unknown time; reporting USR" << endl;
    timeFunction(_USRtime);
    break;
  }
}

double 
OpTimer::tic(unsigned N, unsigned reportInterval) 
{
  _N        = N;
  _interval = reportInterval;
  _i        = 0;
  return (_start = _time());
}    

double 
OpTimer::toc(unsigned i) 
{
  double dt = _time() - _start;
  _i += i;

  if (_verbose && _operation && !(_i % _interval)) {
    _os->setf(ios::fixed);
    int p = _os->precision(3);

    *_os << _TIME_STRINGS[_timeType] << " time elapsed in " << _operation << ": ";
    _printTime(dt);
    if (_N) {
      double fraction = double(_i)/_N;
      *_os << " (" << ROUND(fraction*100) << "% of ";
      _printTime(dt/fraction) << ")";
    }
    *_os << endl;

    _os->precision(p);
    _os->unsetf(ios::fixed);
  }
  
  return dt;
}

double 
OpTimer::_CPUtime()
{
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return (double) usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec / 1.0e6);
}

double 
OpTimer::_SYStime()
{
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return (double) usage.ru_stime.tv_sec + (usage.ru_stime.tv_usec / 1.0e6);
}

double 
OpTimer::_USRtime()
{
  struct timeval T;
  gettimeofday(&T);
  return (double) T.tv_sec + (T.tv_usec / 1.0e6);
}

void 
OpTimer::_newOperation(const char *operation) 
{
  if (_operation) {
    delete [] _operation;
    _operation = 0;
  }
  
  if (operation) {
    _operation = new char [strlen(operation) + 1];
    assert(operation);
    strcpy(_operation, operation);
  }
}

ostream& 
OpTimer::_printTime(double sec) const 
{
  if (_verbose) {
    unsigned hrs = 0;
    
    if (sec >= 3600) {
      hrs = unsigned(sec) / 3600;
      *_os << hrs << ":";
      sec -= hrs * 3600;
    }
    
    if (hrs || (sec >= 60)) {
      unsigned min = unsigned(sec) / 60;
      *_os << min << ":";
      sec -= min * 60;
    }
    
    *_os << sec;
  }

  return *_os;
}

