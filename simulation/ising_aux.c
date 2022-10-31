#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "ising_aux.h"

double correlate(double *x, double *y, unsigned n){
	double sum = 0;
	unsigned i;
	for(i = 0; i < n; i++) sum += x[i]*y[i];
	return sum;
}

void print_vec(double *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\n", vec[i]);
}

void print_cvec(double complex *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\t%+g i\n", creal(vec[i]), cimag(vec[i]));
}

void print_mat(double *m, unsigned long n){
	unsigned i, k;

	for(i = 0; i < n; i++){
		for(k = 0; k < n; k++, m++) printf("%g,\t", *m);
		printf("\n");
	}
}

void print_spins(short *s, unsigned L, unsigned Nt){
	unsigned i, k;

	for(i = 0; i < L; i++){
		for(k = 0; k < Nt; k++) printf("%d\t", s[i*Nt+k]);
		printf("\n");
	}
}
