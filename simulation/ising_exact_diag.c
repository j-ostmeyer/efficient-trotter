#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

#ifndef MKL_LAPACKE
#include <lapacke.h>
#else
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

#include "ising_aux.h"

void mat_mul(double *m, double complex *x, double complex *y, unsigned long n){
	memset(y, 0, n * sizeof(double complex));
	for(unsigned i = 0; i < n; i++)
		for(unsigned k = 0; k < n; k++, m++) y[i] += (*m) * x[k];
}

void mat_mul_tr(double *m, double complex *x, double complex *y, unsigned long n){
	memset(y, 0, n * sizeof(double complex));
	for(unsigned i = 0; i < n; i++)
		for(unsigned k = 0; k < n; k++, m++) y[k] += (*m) * x[i];
}

void set_hamiltonian_heisenberg(double *H, double *J, double *h, unsigned L, unsigned long N){
	for(unsigned i = 0; i < L; i++){
		for(unsigned k = 0; k < N; k++){
			// sigma^x_i sigma^x_{i+1}
			const unsigned l = k ^ (1 << i) ^ (1 << ((i+1)%L));
			H[N*k+l] -= J[0];

			// sigma^y_i sigma^y_{i+1}
			const int sign = ((k >> i) & 1) ^ ((k >> ((i+1)%L)) & 1);
			H[N*k+l] -= sign? J[1]:-J[1];

			// sigma^z_i sigma^z_{i+1}
			H[(N+1)*k] -= sign? -J[2]:J[2];

			// h_i sigma^z_i
			const int signh = (k >> i) & 1;
			H[(N+1)*k] -= signh? -h[i]:h[i];
		}
	}
}

void diagonalise_heisenberg(double *ev, double *hamilton, double *J, double *h, unsigned L, unsigned long N){
	memset(hamilton, 0, N*N * sizeof(double));
	set_hamiltonian_heisenberg(hamilton, J, h, L, N);
	//print_mat(hamilton, N);

	LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, hamilton, N, ev);

	//print_vec(ev, N);
	//print_mat(hamilton, N);
}

void exact_time_evol(double complex *x, double complex *y, double *ev, double *hamilton, double complex t, unsigned long N){
	mat_mul(hamilton, x, y, N);
	for(unsigned i = 0; i < N; i++) y[i] *= cexp(-t * ev[i]);
	mat_mul_tr(hamilton, y, x, N);
}
