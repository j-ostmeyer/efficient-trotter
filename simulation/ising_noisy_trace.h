#ifndef NOISY_TRACE
#define NOISY_TRACE

void trotter_time_step(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int scheme, int first_all_x, unsigned long i);

void trotter_error(double complex *trace, unsigned long index, double complex *x, double complex *y, double complex *z, double *ev, double *hamilton, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x);
void trace_term_trotter(double complex *trace, unsigned long index, double complex *x, double complex *y, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x);

void trace_operator(double complex *trace, unsigned long index, double complex *x, double complex *y, double complex *z, double *ev, double *hamilton, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x, int operator_id);

unsigned long xor_collapse(unsigned long n);
int hadamard_entry(unsigned long i, unsigned long k, unsigned long N);
void hadamard_vector(double complex *x, unsigned long i, unsigned long N);
unsigned long position(unsigned long *indices, unsigned long i, unsigned long N);

double complex *trace_estimator(double *J, double *h, double complex t, unsigned L, unsigned long n_steps, unsigned n, int scheme, int time_series, int full_data, int operator_id, int first_all_x, int seed);
#endif
