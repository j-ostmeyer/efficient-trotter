#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <string.h>

#include "ising_hamiltonian.h"

void apply_exp_hzL(double complex *x, double *J, double *h, double complex t_step, unsigned i, unsigned L, unsigned long N){
	const unsigned shift1 = (i + 1) % L;

	double complex coupling[4] = {cexp(t_step*(J[2]+h[i])), cexp(t_step*(-J[2]+h[i])), cexp(t_step*(J[2]-h[i])), cexp(t_step*(-J[2]-h[i]))};

#pragma omp parallel for
	for(unsigned long k = 0; k < N; k++){
		const unsigned sign = ((k >> i) ^ (k >> shift1)) & 1;
		const unsigned signh = (k >> i) & 1;

		x[k] *= coupling[sign + 2*signh];
	}
}

void apply_exp_hz(double complex *x, double *J, double *h, double complex t_step, unsigned L, unsigned long N){
	for(unsigned i = 0; i < L; i++){
		apply_exp_hzL(x, J, h, t_step, i, L, N);
	}
}

void apply_exp_hxyL(double complex *x, double complex *y, double *J, int pauliY, double complex t_step, unsigned i, unsigned L, unsigned long N, int copy_back){
	const double complex ch = ccosh(J[pauliY] * t_step);
	const double complex sh = csinh(J[pauliY] * t_step);

	const unsigned shift1 = (i + 1) % L;
	const unsigned long pos = 1 << i;
	const unsigned long pos1 = 1 << shift1;

#pragma omp parallel for
	for(unsigned long k = 0; k < N; k++){
		const unsigned long l = k ^ pos ^ pos1;
		const int sign = pauliY? ((((k >> i) ^ (k >> shift1)) & 1)? 1:-1) : 1; // need additional sign for sigma^y, but none for sigma^x

		y[k] = ch * x[k] + sign * sh * x[l];
	}

	if(copy_back) memcpy(x, y, N * sizeof(double complex));
}

void apply_exp_hxy(double complex *x, double complex *y, double *J, int pauliY, double complex t_step, unsigned L, unsigned long N){
	// x: in/out, y: auxiliary only
	for(unsigned i = 0; i < L; i++){
		apply_exp_hxyL(x, y, J, pauliY, t_step, i, L, N, 0);

		double complex *dummy = x;
		x = y;
		y = dummy;
	}

	if(L % 2 == 1) memcpy(y, x, N * sizeof(double complex));
}

void apply_exp_h_pauliL(double complex *x, double complex *y, double *J, double *h, int pauli, double complex t_step, unsigned i, unsigned L, unsigned long N){
	if(pauli == 3) apply_exp_hzL(x, J, h, t_step, i, L, N);
	else apply_exp_hxyL(x, y, J, pauli-1, t_step, i, L, N, 1);
}

void apply_exp_h_pauli(double complex *x, double complex *y, double *J, double *h, int pauli, double complex t_step, unsigned L, unsigned long N){
	if(pauli == 3) apply_exp_hz(x, J, h, t_step, L, N);
	else apply_exp_hxy(x, y, J, pauli-1, t_step, L, N);
}

void apply_exp_h(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, int forward){
	if(first_all_x){
		if(forward)
			for(unsigned pauli = 1; pauli <= 3; pauli++) apply_exp_h_pauli(x, y, J, h, pauli, t_step, L, N);
		else
			for(unsigned pauli = 3; pauli >= 1; pauli--) apply_exp_h_pauli(x, y, J, h, pauli, t_step, L, N);
	}else{
		if(forward)
			for(unsigned i = 0; i < L; i++)
				for(unsigned pauli = 1; pauli <= 3; pauli++) apply_exp_h_pauliL(x, y, J, h, pauli, t_step, i, L, N);
		else
			for(int i = L-1; i >=0; i--)
				for(unsigned pauli = 3; pauli >= 1; pauli--) apply_exp_h_pauliL(x, y, J, h, pauli, t_step, i, L, N);
	}
}

void apply_hamiltonian(double complex *x, double complex *y, double *J, double *h, unsigned L, unsigned long N){
	memset(y, 0, N * sizeof(double complex));

	for(unsigned i = 0; i < L; i++){
		const unsigned shift1 = (i + 1) % L;
		const unsigned long pos = 1 << i;
		const unsigned long pos1 = 1 << shift1;

#pragma omp parallel for
		for(unsigned long k = 0; k < N; k++){
			// sigma^x_i sigma^x_{i+1} + sigma^y_i sigma^y_{i+1}
			const unsigned long l = k ^ pos ^ pos1;
			const unsigned sign = ((k >> i) ^ (k >> shift1)) & 1;
			const double complex termXY = (J[0] + (sign? J[1]:-J[1])) * x[l];

			// sigma^z_i sigma^z_{i+1} + h_i sigma^z_i
			const unsigned signh = (k >> i) & 1;
			const double couplingZ = (sign? -J[2]:J[2]) + (signh? -h[i]:h[i]);
			y[k] -= couplingZ * x[k] + termXY;
		}
	}
}
