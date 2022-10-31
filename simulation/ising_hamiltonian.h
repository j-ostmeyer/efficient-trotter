#ifndef ISING_HAMILTONIAN
#define ISING_HAMILTONIAN

void apply_exp_hzL(double complex *x, double *J, double *h, double complex t_step, unsigned i, unsigned L, unsigned long N);
void apply_exp_hz(double complex *x, double *J, double *h, double complex t_step, unsigned L, unsigned long N);

void apply_exp_hxyL(double complex *x, double complex *y, double *J, int pauliY, double complex t_step, unsigned i, unsigned L, unsigned long N, int copy_back);
void apply_exp_hxy(double complex *x, double complex *y, double *J, int pauliY, double complex t_step, unsigned L, unsigned long N);

void apply_exp_h_pauliL(double complex *x, double complex *y, double *J, double *h, int pauli, double complex t_step, unsigned i, unsigned L, unsigned long N);
void apply_exp_h_pauli(double complex *x, double complex *y, double *J, double *h, int pauli, double complex t_step, unsigned L, unsigned long N);
void apply_exp_h(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, int forward);

void apply_hamiltonian(double complex *x, double complex *y, double *J, double *h, unsigned L, unsigned long N);

void magnetisation(double complex *x, double *m, unsigned L, unsigned long N);
#endif
