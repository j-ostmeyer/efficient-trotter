#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <time.h>

#include "mt19937-64.h"
#include "ising_aux.h"
#include "ising_hamiltonian.h"
#include "ising_trotter.h"
#include "ising_exact_diag.h"
#include "ising_noisy_trace.h"

void trotter_time_step(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int scheme, int first_all_x, unsigned long i){
	switch(scheme){ // digits: order, cycles, id
		case 211:
			leap_frog(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 221:
			omelyan2(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 431:
			forest_ruth(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 441:
			fr_type(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 442:
			smallB(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 443:
			non_unitary1(x, y, J, h, t_step, L, N, first_all_x, i);
			break;
		case 451:
			st4(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 452:
			alg30(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 453:
			opt_ord4(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 454:
			non_unitary2(x, y, J, h, t_step, L, N, first_all_x, i);
			break;
		case 455:
			omelyan_st4(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 456:
			non_unitary3(x, y, J, h, t_step, L, N, first_all_x, i);
			break;
		case 457:
			non_unitary_blanes(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 458:
			non_unitary_const(x, y, J, h, t_step, L, N, first_all_x, i);
			break;
		case 459:
			close_opt4_unif(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 461:
			blanes4(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 471: // technically this one is only order 4 (it's order 6 when some of the commutators vanish)
			order6(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 671:
			yoshida6a(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 691:
			fr_squared(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 6101:
			blanes6(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 6151:
			fr_st6(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 6251:
			st6(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 8171:
			morales8_8(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 8271:
			fr_cubed(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 8451:
			fr_sq_st8(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 8501:
			blanes6_st8(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 8502:
			blanes6_yoshi8(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 81251:
			st8(x, y, J, h, t_step, L, N, first_all_x);
			break;
		case 41: // not polynomial (order 0)!
			taylor4(x, y, J, h, t_step, L, N);
			break;
		case 171: // not polynomial (order 0)!
			taylor17(x, y, J, h, t_step, L, N);
			break;
		case 521: // not polynomial (order 0)!
			taylor52(x, y, J, h, t_step, L, N);
			break;
		case 522: // not polynomial (order 0)!
			taylor52_prod(x, y, J, h, t_step, L, N);
			break;
		case 881: // not polynomial (order 0)!
			taylor88(x, y, J, h, t_step, L, N);
			break;
		case 882: // not polynomial (order 0)!
			taylor88_prod(x, y, J, h, t_step, L, N);
			break;
		case 3041: // not polynomial (order 0)!
			taylor304(x, y, J, h, t_step, L, N);
			break;
		case 3042: // not polynomial (order 0)!
			taylor304_prod(x, y, J, h, t_step, L, N);
			break;
		case 1521: // not polynomial (order 0)!
			chebyshev152(x, y, J, h, t_step, L, N);
			break;
		case 1522: // not polynomial (order 0)!
			chebyshev152_prod(x, y, J, h, t_step, L, N);
			break;
		default:
			printf("Decomposition scheme %d not supported!\n", scheme);
			exit(0);
	}
}

void trotter_error(double complex *trace, unsigned long index, double complex *x, double complex *y, double complex *z, double *ev, double *hamilton, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x){
#ifndef SKIP_EXACT_DIAG
	const double complex t_step = t / n_steps;

	hadamard_vector(x, index, N);

	for(unsigned long i = 0; i < n_steps; i++){
		trotter_time_step(x, y, J, h, t_step, L, N, scheme, first_all_x, i);

		if(time_series || i == n_steps-1){
			hadamard_vector(z, index, N);
			exact_time_evol(z, y, ev, hamilton, (i+1)*t_step, N);

			double complex trX = 0, trZ = 0;
			double diff = 0;
			for(unsigned long k = 0; k < N; k++){
				trX += hadamard_entry(index, k, N) * x[k];
				trZ += hadamard_entry(index, k, N) * z[k];
				diff += pow(cabs(z[k] - x[k]), 2);
			}
			if(full_data)
				trace[time_series? i:0] += diff;
			else{
				trace[time_series? i:0] += trX - trZ;
				trace[time_series? n_steps+i:1] += diff;
			}
		}
	}
#endif
}

void trace_term_trotter(double complex *trace, unsigned long index, double complex *x, double complex *y, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x){
	const double complex t_step = t / n_steps;

	hadamard_vector(x, index, N);

	for(unsigned long i = 0; i < n_steps; i++){
		trotter_time_step(x, y, J, h, t_step, L, N, scheme, first_all_x, i);

		if(time_series || i == n_steps-1){
			double complex tr = 0;
#pragma omp parallel for reduction(+:tr)
			for(unsigned long k = 0; k < N; k++) tr += hadamard_entry(index, k, N) * x[k];
			trace[time_series? i:0] += tr;
			if(!full_data)
				trace[time_series? n_steps+i:1] += creal(tr)*creal(tr) + I*cimag(tr)*cimag(tr); // std dev for Re and Im separately
		}
	}
}

void trace_operator(double complex *trace, unsigned long index, double complex *x, double complex *y, double complex *z, double *ev, double *hamilton, double *J, double *h, double complex t, unsigned L, unsigned long N, unsigned long n_steps, int scheme, int time_series, int full_data, int first_all_x, int operator_id){
		switch(operator_id){
			case -1: // compare time evolution to exact result
				trotter_error(trace, index, x, y, z, ev, hamilton, J, h, t, L, N, n_steps, scheme, time_series, full_data, first_all_x);
				break;
			case 0: // compute spectral form factor (SFF)
				trace_term_trotter(trace, index, x, y, J, h, t, L, N, n_steps, scheme, time_series, full_data, first_all_x);
				break;
			case 1: // implement your favourite operator here!
				break;
			default:
				printf("Operator with id %d not known!\n", operator_id);
				exit(0);
		}
}

unsigned long xor_collapse(unsigned long n){
	//returns 1 if n contains an even amount of 1s, 0 else
	n ^= n >> 32;
	n ^= n >> 16;
	n ^= n >> 8;
	n ^= n >> 4;
	n ^= n >> 2;
	n ^= n >> 1;

	return n & 1;
}

int hadamard_entry(unsigned long i, unsigned long k, unsigned long N){
	const unsigned long sign = xor_collapse(i & k);
	return sign? -1:1;
}

void hadamard_vector(double complex *x, unsigned long i, unsigned long N){
#pragma omp parallel for
	for(unsigned long k = 0; k < N; k++){
		x[k] = hadamard_entry(i, k, N);
	}
}

unsigned long position(unsigned long *indices, unsigned long i, unsigned long N){
	if(indices){
		const unsigned long pos = i + genrand64_real2()*(N-i);
		const unsigned long index = indices[pos];
		indices[pos] = indices[i];
		indices[i] = index;
		return index;
	}else return genrand64_real2() * N;
}

double complex *trace_estimator(double *J, double *h, double complex t, unsigned L, unsigned long n_steps, unsigned n, int scheme, int time_series, int full_data, int operator_id, int first_all_x, int seed){
	const unsigned long N = 1L << L;
	const unsigned long out_dim = (time_series? n_steps:1) * (full_data? n:2);
	double complex *x = malloc(N * sizeof(double complex));
	double complex *y = malloc(N * sizeof(double complex));
	double complex *trace = calloc(out_dim, sizeof(double complex));
	unsigned long *indices = NULL;

	double *ev = operator_id<0? malloc(N * sizeof(double)) : NULL;
	double *hamilton = operator_id<0? malloc(N*N * sizeof(double)) : NULL;
	double complex *z = operator_id<0? malloc(N * sizeof(double complex)) : NULL;

#ifndef SKIP_EXACT_DIAG
	if(operator_id<0) diagonalise_heisenberg(ev, hamilton, J, h, L, N);
#endif

	//if(scheme == 41 || scheme == 171 || scheme == 521 || scheme == 881 || scheme == 3041)
		y = realloc(y, 2*N * sizeof(double complex));

	if(scheme == 1521)
		y = realloc(y, 3*N * sizeof(double complex));

	init_genrand64(time(NULL)+seed);

	if(n > sqrt(N)){ // only sample without repetition when it's worth it
		if(n > N) n = N;
		indices = malloc(N * sizeof(unsigned long));
#pragma omp parallel for
		for(unsigned long k = 0; k < N; k++) indices[k] = k;
	}

	for(unsigned long i = 0; i < n; i++){
		const unsigned long index = position(indices, i, N);
		const unsigned long offset = i * full_data * (time_series? n_steps:1);

		trace_operator(trace+offset, index, x, y, z, ev, hamilton, J, h, t, L, N, n_steps, scheme, time_series, full_data, first_all_x, operator_id);
	}

	if(!full_data){
		const unsigned long n_times = time_series? n_steps:1;
		for(unsigned long i = 0; i < n_times; i++){
			if(operator_id < 0){
				trace[n_times + i] /= n;
				trace[n_times + i] = sqrt(creal(trace[n_times + i]));
			}else{
				trace[n_times + i] -= (pow(creal(trace[i]),2) + I*pow(cimag(trace[i]),2)) / n;
				trace[n_times + i] /= n-1;
				trace[n_times + i] = sqrt(creal(trace[n_times + i])) + I*sqrt(cimag(trace[n_times + i]));
			}
			trace[i] /= n;
		}
	}

	free(x);
	free(y);
	if(indices) free(indices);

	if(ev) free(ev);
	if(hamilton) free(hamilton);
	if(z) free(z);

	return trace;
}
