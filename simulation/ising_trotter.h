#ifndef ISING_TROTTER
#define ISING_TROTTER

void print_coeffs(int cycles, double complex *a, double complex *b);

void decompose_trotter(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int cycles, double complex *a, double complex *b, int first_all_x);
void decompose(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int cycles, double complex *a, double complex *b, int first_all_x);

void leap_frog(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void omelyan2(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void forest_ruth(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void fr_type(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void smallB(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void non_unitary1(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i);

void st4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void alg30(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void opt_ord4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void non_unitary2(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i);
void omelyan_st4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void non_unitary3(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i);
void non_unitary_blanes(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);
void non_unitary_const(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i);
void close_opt4_unif(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void blanes4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void order6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void fr_squared(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void blanes6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void fr_st6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void st6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void fr_cubed(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void fr_sq_st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void blanes6_st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x);

void taylor17(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N);
#endif
