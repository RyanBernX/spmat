CC=gcc
CFLAG=-I../include -I$(MKLROOT)/include -fopenmp
CDEBUG=-g -Wall
COPTFLAG=-O2
LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
EOBJ=main.o power_impls.o spmat.o
LOBJ=power_impls.o spmat.o


%.o: %.c
	$(CC) $(CDEBUG) $(COPTFLAG) -c $< $(CFLAG)

.PHONY: lib exec clean

lib: $(LOBJ)
	ar rv libeigpow.a $(LOBJ)
	ranlib libeigpow.a

exec: $(EOBJ)
	$(CC) $(CDEBUG) $(COPTFLAG) -o main.x $(EOBJ) $(LIBS) $(CFLAG)

clean:
	rm -f *.o *.x *.a

