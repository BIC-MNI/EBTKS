CFLAGS = -O2

OBJS = dsysv.o dsytf2.o dsytrf.o dsytrs.o dlasyf.o \
     idamax.o dcopy.o dgemm.o dgemv.o dger.o dscal.o dswap.o dsyr.o \
     ieeeck.o ilaenv.o lsame.o xerbla.o \
     s_cmp.o s_copy.o

all: libdsysv.a

libdsysv.a: $(OBJS)
	ar r libdsysv.a $(OBJS)
	ranlib libdsysv.a

$(OBJS): f2c.h blaswrap.h

