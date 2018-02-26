# Makefile for aeden. T.Pope, 2018
FC          = gfortran 
FLAGS_DEBUG = -g
FLAGS_OPT   = -O3
FLAGS_FORM  = -ffixed-line-length-none 
LAPACK      = /usr/lib/liblapack.so.3
BLAS        = /usr/lib/libblas.so.3 -fexternal-blas
LIBRARY     = $(BLAS) $(LAPACK) 

FLAGS       = $(FLAGS_DEBUG) $(FLAGS_OPT) $(FLAGS_FORM) 

.f.o:
	$(FC) -c $(FLAGS) $<

OBJS = rundata.o functions.o basis.o scf.o io.o aeden.o

ai:	$(OBJS)
	$(FC) -o aeden.x $(OBJS) $(FLAGS) $(LIBRARY)

clean:
	rm *.o; rm *.mod
