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
$RCSfile: popen.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#include "popen.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream.h>
#include <osfcn.h>
#include <signal.h>

#if SVR4 || __STDC__ || sun
extern "C" pid_t waitpid(pid_t, int*, int);

#define App(x)

inline pid_t Wait(pid_t p, int *x=0) { return waitpid(p, x, 0); }

#else
extern "C" pid_t wait(int*);

#define App(p)  plist.append(new Pid(p))

// unsecure wait, solved by linked list administration of child pids
#include <tlist.h>
class Pid{
public:
	Pid(pid_t p) : pd(p), rg(1), extcod(-1){}
	pid_t pid() { return pd; }
	int   exitcode() { return extcod; }
	void notrunning(int x) { rg = 0; extcod = x ;}
	int running() { return rg;}
private:
	pid_t pd;
	int extcod;
	int rg;
};

static tlist<Pid> plist;

// Wait() does not check for stopped processes.
pid_t Wait(pid_t p, int *fx=0)
{
	// look up p, remember place:
	int lx;
	Pid *pq;
	int cq;
	for (cq = 0; (pq = plist.get(cq)) && (pq->pid() != p); ++cq)
		; // void

	if (!pq){
		errno = ESRCH;
		return -1;
	}

	// stopped before ?
	if (!pq->running()){
		if (fx) *fx = pq->exitcode();
		plist.remove(cq);
		delete pq;
		return p;
	}

	// look up the proces by waiting until found
	while (1){
		pid_t rp = wait(&lx);

		if (rp == -1){
			errno = ECHILD;
			return -1;
		}

		if (p == rp){
			// got a match, remove from plist
			if (fx) *fx = lx;
			plist.remove(cq);
			delete pq;
			return p;
		}

		// not the proces wanted, continue, but mark as
		// being removed;  look up the proces
		Pid *rq;
		int rc;	
		for (rc = 0; (rq = plist.get(rc)) && (rq->pid() != rp);
							++rc)
			; // void

		if (!rq){
			errno = ESRCH;
			return -1;
		}

		rq->notrunning(lx);
	}
}

#endif

int
psystem(const char *command)
{
  long pid;
  int ex;

	switch (pid = ::fork()){
	case 0:
		// child
		::putenv("IFS=");
		::setuid(getuid());
		::setgid(getgid());
		::execl("/bin/sh", "sh", "-c", command, 0);
		::kill(getpid(), SIGUSR1);	//inform parent
	case -1:
		return -1;
	default:
		// parent
		App(pid);
		if (Wait(pid, &ex) != pid)
			return -1;

		if ((ex&255) == SIGUSR1){
			errno = EACCES;
			return -1;
		}
		return ex;
	}
}

int
popenbuf::overflow(int c)
{
	int df = pptr()-pbase();
	int n  = c;

	if (df > 0 && write(pfdout, pbuf, df) != df)
		n = -1;

	if (c != EOF && n != -1){
		char _c = c;
		if (write(pfdout, &_c, 1) != 1)
			n = -1;
	}

	setp(base(), ebuf());

	return n == -1 ? -1 : zapeof(c);
}

int
popenbuf::underflow()
{
	// flush output
	if (pfdout >= 0) overflow();

	int n = read(pfdin, gbuf, bufsize);

	if (n < 0){
		setg(0, 0, 0);
		return EOF;
	}

	if (n == 0) 
		setg(gbuf, gbuf, gbuf+1);
	else
		setg(gbuf, gbuf, gbuf+n);

	return n == 0 ? EOF : zapeof(*gbuf);
}

//static char *getb(int n){ return new char[n];}

popenbuf::popenbuf(const char *command, ios::open_mode om) :
  pbuf(new char[bufsize*2]),
  gbuf(pbuf+bufsize), 
  running(0), 
  opm(om), 
  pfdin(-1),
  pfdout(-1), 
  pid(-1)
{
  setbuf(pbuf, bufsize),

	pbuf = base();

	assert(pbuf);

	setp(base(), ebuf());
	setg(0, 0, 0);

	switch(opm){
	case ios::in:
	case ios::out:
		popen1(command);
		return;
	case ios::in|ios::out:
		popen2(command);
		return;
	default:
		return;		// throw.....
	};
}

void popenbuf::popen1(const char *command)
{
	// one way pipe

	int pip[2];

	// set up pipe for child
	if (::pipe(pip) < 0){
		return	;	// throw
	}

	switch (pid=::fork()){
	case 0:
		// child
		if (opm == ios::in){
			::close(1);
			::dup(pip[1]);
		} else {
			::close(0);
			::dup(pip[0]);
		}
		::close(pip[0]);
		::close(pip[1]);
		::putenv("IFS=");
		::setuid(getuid());
		::setgid(getgid());
		::execl("/bin/sh", "sh", "-c", command, 0);
		::kill(getpid(), SIGUSR1);	//inform parent
	case -1:
		// failed
		return;		// throw
	default:
		// parent
		running = 1;
		App(pid);
		if (opm == ios::in){
			::close(pip[1]);
			pfdin = pip[0];
		} else {
			::close(pip[0]);
			pfdout = pip[1];
		}
	}
}

void
popenbuf::popen2(const char *command)
{
	// set up 2 pipes

	int pws[2];
	int prs[2];

	if (::pipe(pws) < 0){
		return	;	// throw
	}

	if (::pipe(prs) < 0){
		::close(pws[0]);
		::close(pws[1]);
		return	;	// throw
	}

	switch (pid=::fork()){
	case 0:
	  cout << "0 " << flush;
		// child
		::close(0);      ::dup(pws[0]);
		::close(pws[0]); ::close(pws[1]);

		::close(1);      ::dup(prs[1]);
		::close(prs[0]); ::close(prs[1]);

		::putenv("IFS=");
		::setuid(getuid());
		::setgid(getgid());
		::execl("/bin/sh", "sh", "-c", command, 0);
		::kill(getpid(), SIGUSR1);	//inform parent
	case -1:
	  cout << "-1 " << flush;
		// failed
		return;		// throw
	default:
	  cout << "default " << flush;
		// parent
		running = 1;
		App(pid);
		::close(pws[0]);
		::close(prs[1]);
		pfdout = pws[1];
		pfdin  = prs[0];
	}
}

popenbuf::~popenbuf()
{
	exit();
	delete pbuf;
	pbuf = gbuf = 0;
}

int
popenbuf::exit()
{
	int ex;

	if (!running){
		errno = ESRCH;
		return -1;
	}

	if ((opm&ios::out) && pfdout >= 0)
		overflow();

	if (pfdout >= 0) ::close(pfdout);	
	if (pfdin  >= 0) ::close(pfdin);	

	if (::Wait(pid, &ex) != pid)
		return -1;

	running = 0;
	pid = -1;
	opm=ios::open_mode(0);
	pfdin = pfdout = -1;

	if ((ex&255) == SIGUSR1){
		errno = EACCES;
		return -1;
	}

	return ex;
}

void
popenbuf::closewriteside()
{
	if (!running)
		return;

	if (opm != (ios::in|ios::out))
		return;

	if (pfdout >= 0) overflow();
	::close(pfdout);
	pfdout = -1;
}

pstreambase::pstreambase(const char *c, ios::open_mode om) : buf(c, om){ init(&buf);}
int pstreambase::exit() { return buf.exit(); }
void pstreambase::closewriteside() { buf.closewriteside(); }

opopen::opopen(const char *c) : pstreambase(c, ios::out) { }
int opopen::exit() { return pstreambase::exit(); }

ipopen::ipopen(const char *c) : pstreambase(c, ios::in) { }
int ipopen::exit() { return pstreambase::exit(); }

iopopen::iopopen(const char *c) : pstreambase(c, ios::open_mode(ios::in|ios::out)){}
int iopopen::exit() { return pstreambase::exit(); }
void iopopen::closewriteside() { pstreambase::closewriteside(); }
