MYFLAGS=-Wall -pedantic -O2 -flto -fPIC #-g -Og -fsanitize=address -static-libasan
MYLIBS=-lm
CC = gcc
MYFILES=aux hamiltonian trotter noisy_trace
RNG=mt19937-64

MYFLAGS += -Wno-unused-result

ifneq ($(filter $(HOSTNAME),d3-login02 d3-login01 viz01.pri.barkla.alces.network),)
#  #We are on DIaL
  #CC = icc -std=c11
  #MYFLAGS += -qopenmp -pthread -no-multibyte-chars
  MYFLAGS += -fopenmp
  #MYFLAGS += -I/opt/gridware/depots/e2b91392/el7/pkg/apps/R/4.1.0/gcc-5.5.0+lapack-3.5.0+blas-3.6.0/lib64/R/include
  MYFLAGS += -I/opt/gridware/depots/e2b91392/el7/pkg/apps/R/3.6.3/gcc-5.5.0+lapack-3.5.0+blas-3.6.0/lib64/R/include
  MYFLAGS += -DSKIP_EXACT_DIAG
else
  MYFILES += exact_diag
  MYFLAGS += -std=gnu11 -march=native
  MYFLAGS += -I/usr/share/R/include
  MYFLAGS += -Wno-unknown-pragmas
  MYLIBS += -llapacke -llapack -lblas
endif

all: ising_dos.so

# Take file identifiers and change empty prefix by ising_ % empty suffix by .o
# Compile all dependencies $^ to target $@

ising_dos.so: $(MYFILES:%=ising_%.o) $(RNG:%=%.o) ising_R.o
	$(CC) $(MYFLAGS) -shared -o $@ $^ $(MYLIBS)

# For every object .o check for .c and .h file with same name
# Compile first dependency $< (the .c file)
%.o: %.c %.h Makefile
	$(CC) $(MYFLAGS) -c $<

clean:
	rm -f *.o *.so ising_dos ising_noisy_trace exact

distclean:
	rm -f *.o *.so *.exe *~ tags ising_dos ising_noisy_trace exact
