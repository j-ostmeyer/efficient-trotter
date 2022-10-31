#ifndef EXACT_DIAG
#define EXACT_DIAG

void mat_mul(double *m, double complex *x, double complex *y, unsigned long n);
void mat_mul_tr(double *m, double complex *x, double complex *y, unsigned long n);

void set_hamiltonian_heisenberg(double *H, double *J, double *h, unsigned L, unsigned long N);
void diagonalise_heisenberg(double *ev, double *hamilton, double *J, double *h, unsigned L, unsigned long N);

void exact_time_evol(double complex *x, double complex *y, double *ev, double *hamilton, double complex t, unsigned long N);

#endif
