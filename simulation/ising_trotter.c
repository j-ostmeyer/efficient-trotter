#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <time.h>

#include "ising_hamiltonian.h"
#include "ising_trotter.h"

void print_coeffs(int cycles, double complex *a, double complex *b){
	int i;
	cycles = abs(cycles);
	printf("\\begin{align}\n");
	printf("\\begin{split}\n");
	for(i = 0; i < (cycles-1)/2; i++){
		if(cimag(a[i]) == 0 && cimag(b[i]) == 0)
			printf("a_{%d} &= %.16g\\,, \\\\ b_{%d} &= %.16g\\,,\\\\\n", i+1, creal(b[i]), i+1, creal(a[i]));
		else
			printf("a_{%d} &= %.16g %+.16g \\im\\,, \\\\ b_{%d} &= %.16g %+.16g \\im\\,,\\\\\n", i+1, creal(b[i]), cimag(b[i]), i+1, creal(a[i]), cimag(a[i]));
	}

	if(cycles%2 == 0){
		if(cimag(b[i]) == 0)
			printf("a_{%d} &= %.16g\\,, \\\\ b_{%d} &= \\frac 12 - \\sum_{i=1}^{%d}b_i\\,,\\\\\n", i+1, creal(b[i]), i+1, i);
		else
			printf("a_{%d} &= %.16g %+.16g \\im\\,, \\\\ b_{%d} &= \\frac 12 - \\sum_{i=1}^{%d}b_i\\,,\\\\\n", i+1, creal(b[i]), cimag(b[i]), i+1, i);

		i++;
		printf("a_{%d} &= 1 - 2\\sum_{i=1}^{%d}a_i\\,.\n", i+1, i);
	}else printf("a_{%d} &= \\frac 12 - \\sum_{i=1}^{%d}a_i\\,, \\\\ b_{%d} &= 1 - 2\\sum_{i=1}^{%d}b_i\\,.\n", i+1, i, i+1, i);

	printf("\\end{split}\\label{eq:%d-}\n", cycles);
	printf("\\end{align}\n\n");
}

void decompose_trotter(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int cycles, double complex *a, double complex *b, int first_all_x){
	double complex c = b[0];

	for(int i = 0; i < abs(cycles)/2; i++){
		apply_exp_h(x, y, J, h, c * t_step, L, N, first_all_x, 1);
		c = a[i] - c;
		apply_exp_h(x, y, J, h, c * t_step, L, N, first_all_x, 0);
		c = b[i+1] - c;
	}

	for(int i = (abs(cycles)-1)/2; i >= 0; i--){
		apply_exp_h(x, y, J, h, c * t_step, L, N, first_all_x, 1);
		c = a[i] - c;
		apply_exp_h(x, y, J, h, c * t_step, L, N, first_all_x, 0);
		c = (cycles > 0? b[i]: conj(b[i])) - c;
	}
}

void decompose(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int cycles, double complex *a, double complex *b, int first_all_x){
	decompose_trotter(x, y, J, h, t_step, L, N, cycles, a, b, first_all_x);

	//print_coeffs(cycles, a, b);
}

void leap_frog(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	double complex a[1] = {1};
	double complex b[1] = {.5};

	decompose(x, y, J, h, t_step, L, N, 1, a, b, first_all_x);
}

void omelyan2(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double lambda = 0.1931833275037836;
	double complex a[1] = {.5};
	double complex b[2] = {lambda, 1-2*lambda};

	decompose(x, y, J, h, t_step, L, N, 2, a, b, first_all_x);
}

void forest_ruth(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double lambda = 0.6756035959798288, theta = 2*lambda;
	double complex a[2] = {theta, 1-2*theta};
	double complex b[2] = {lambda, .5-lambda};

	decompose(x, y, J, h, t_step, L, N, 3, a, b, first_all_x);
}

void fr_type(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double xi = 0.1720865590295143, lambda = -0.09156203075515678, chi = -0.1616217622107222;
	double complex a[2] = {.5-lambda, lambda};
	double complex b[3] = {xi, chi, 1-2*(xi+chi)};

	decompose(x, y, J, h, t_step, L, N, 4, a, b, first_all_x);
}

void smallB(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double xi = 0.5316386245813512, lambda = 0.5437514219173741, chi = -0.3086019704406066;
	double complex a[2] = {.5-lambda, lambda};
	double complex b[3] = {xi, chi, 1-2*(xi+chi)};

	decompose(x, y, J, h, t_step, L, N, 4, a, b, first_all_x);
}

void non_unitary1(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i){
	const double complex xi = 0.09957801119428373 + 0.02359386141367452*I;
	const double complex lambda = 0.2403781402426499 - 0.08909472525370253*I;
	const double complex chi = 0.2520542187700347 + 0.09826170579213035*I;

	double complex a[2] = {.5-lambda, lambda};
	double complex b[3] = {xi, chi, 1-2*(xi+chi)};

	if(i%2){
		for(int k = 0; k < 2; k++){
			a[k] = conj(a[k]);
			b[k] = conj(b[k]);
		}
		b[2] = conj(b[2]);
	}

	decompose(x, y, J, h, t_step, L, N, 4, a, b, first_all_x);
}

void st4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s = 0.4144907717943757;
	double complex a[3] = {s, s, 1-4*s};
	double complex b[3] = {.5*s, s, .5-1.5*s};

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void alg30(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double rho = 0.08398315262876693, theta = 0.6822365335719091;
	const double varT = 0.2539785108410595, lambda = -0.03230286765269967;

	double complex a[3] = {varT, lambda, 1-2*(varT+lambda)};
	double complex b[3] = {rho, theta, .5-rho-theta};

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void opt_ord4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double rho = 0.09257547473195787, theta = 0.4627160310210738;
	const double varT = 0.2540996315529392, lambda = -0.1676517240119692;

	double complex a[3] = {varT, lambda, 1-2*(varT+lambda)};
	double complex b[3] = {rho, theta, .5-rho-theta};

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void non_unitary2(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i){
	const double complex rho = 0.07613272445178274 - 0.03518797331257356*I;
	const double complex theta = 0.2017183745725757 + 0.02597491015915232*I;
	const double complex varT = 0.1658339349217486 - 0.07090293766092534*I;
	const double complex lambda = 0.2137425142256234 + 0.1386193640914034*I;

	double complex a[3] = {varT, lambda, 1-2*(varT+lambda)};
	double complex b[3] = {rho, theta, .5-rho-theta};

	if(i%2) for(int k = 0; k < 3; k++){
		a[k] = conj(a[k]);
		b[k] = conj(b[k]);
	}

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void omelyan_st4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double xi = 0.3221375960817984;
	const double lambda = 0.5413165481700430;

	double complex a[3] = {xi, lambda, 1-2*(xi+lambda)};
	double complex b[3] = {.5*xi, .5*(xi+lambda), .5-xi-.5*lambda};

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void non_unitary3(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i){
	//const double complex xi = 0.1558830917220528 + 0.08079680289595839*I;
	//const double complex lambda = 0.2148919565328240 - 0.1432835769224897*I;
	const double complex xi = 0.1817097484549141-0.1080260385004353*I;
	const double complex lambda = 0.2033014198470074+0.1428634440888616*I;
	//const double complex xi = 0.1927546141028121 + 0.0552778546008495*I, lambda = xi;

	double complex a[3] = {xi, lambda, 1-2*(xi+lambda)};
	double complex b[3] = {.5*xi, .5*(xi+lambda), .5-xi-.5*lambda};

	if(i%2) for(int k = 0; k < 3; k++){
		a[k] = conj(a[k]);
		b[k] = conj(b[k]);
	}

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void non_unitary_blanes(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	double complex a[3] = {1./8, 0.23670501659941197298, 0.27658996680117605403};
	double complex b[3] = {0.03881396214419327198 - 0.045572109263923104872*I, 0.19047619047619047619 + 0.115462072300408741306*I, 0.27070984737961625182 - 0.148322245509626403888*I};

	decompose(x, y, J, h, t_step, L, N, -5, a, b, first_all_x);
}

void non_unitary_const(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x, unsigned long i){
	const double complex rho = 0.1 + 0.02523113193557069*I;
	const double complex theta = 0.2 - 0.04082482904638631*I;
	const double complex varT = 0.2 + 0.05046226387114138*I;
	const double complex lambda = 0.2 - 0.132111921963914*I;

	double complex a[3] = {varT, lambda, 1-2*(varT+lambda)};
	double complex b[3] = {rho, theta, .5-rho-theta};

	if(i%2) for(int k = 0; k < 3; k++){
		a[k] = conj(a[k]);
		b[k] = conj(b[k]);
	}

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void close_opt4_unif(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double rho = 0.0935624228832911, theta = 0.4905993272009264;
	const double varT = 0.2570861302065004, lambda = -0.113194296843695;

	double complex a[3] = {varT, lambda, 1-2*(varT+lambda)};
	double complex b[3] = {rho, theta, .5-rho-theta};

	decompose(x, y, J, h, t_step, L, N, 5, a, b, first_all_x);
}

void blanes4(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	double complex a[3] = {0.209515106613362, -0.143851773179818, .5-0.209515106613362+0.143851773179818};
	double complex b[4] = {0.0792036964311957, 0.353172906049774, -0.0420650803577195, 1-2*(0.0792036964311957+0.353172906049774-0.0420650803577195)};

	decompose(x, y, J, h, t_step, L, N, 6, a, b, first_all_x);
}

void order6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	double complex a[4] = {0.2465881872786138, 0.6047073875057809, -0.4009869039788007, 0.9938265838881204e-1};
	double complex b[4] = {0.8333333333333333e-1, 0.3977675859548440, -0.3933369314462574e-1, 0.5823277385644840e-1};

	decompose(x, y, J, h, t_step, L, N, 7, a, b, first_all_x);
}

void fr_squared(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s4 = 1.174671758089363;

	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*(1-2*s4), L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
}

void blanes6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	double complex a[5] = {0.148816447901042, -0.132385865767784, 0.067307604692185, 0.432666402578175, .5-0.148816447901042+0.132385865767784-0.067307604692185-0.432666402578175};
	double complex b[6] = {0.0502627644003922,  0.413514300428344, 0.0450798897943977, -0.188054853819569, 0.541960678450780, 1-2*(0.0502627644003922+0.413514300428344+0.0450798897943977-0.188054853819569+0.541960678450780)};

	decompose(x, y, J, h, t_step, L, N, 10, a, b, first_all_x);
}

void fr_st6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s4 = 0.3730658277332728;

	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*(1-4*s4), L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
	forest_ruth(x, y, J, h, t_step*s4, L, N, first_all_x);
}

void st6(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s2 = 0.4144907717943757, s4 = 0.3730658277332728;
	double complex a[13] = {s2*s4, s2*s4, (1-4*s2)*s4, s2*s4, s2*s4, s2*s4, s2*s4, (1-4*s2)*s4, s2*s4, s2*s4, s2*(1-4*s4), s2*(1-4*s4), (1-4*s2)*(1-4*s4)};
	double complex b[13] = {.5*s2*s4, s2*s4, (.5-1.5*s2)*s4, (.5-1.5*s2)*s4, s2*s4, s2*s4, s2*s4, (.5-1.5*s2)*s4, (.5-1.5*s2)*s4, s2*s4, .5*s2*(1-3*s4), s2*(1-4*s4), (.5-1.5*s2)*(1-4*s4)};

	decompose(x, y, J, h, t_step, L, N, 25, a, b, first_all_x);
}

void fr_cubed(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s6 = 1.116182939325386;

	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*(1-2*s6), L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
}

void fr_sq_st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s6 = 0.3595846493499923;

	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*(1-4*s6), L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
	fr_squared(x, y, J, h, t_step*s6, L, N, first_all_x);
}

void blanes6_st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s6 = 0.3595846493499923;

	blanes6(x, y, J, h, t_step*s6, L, N, first_all_x);
	blanes6(x, y, J, h, t_step*s6, L, N, first_all_x);
	blanes6(x, y, J, h, t_step*(1-4*s6), L, N, first_all_x);
	blanes6(x, y, J, h, t_step*s6, L, N, first_all_x);
	blanes6(x, y, J, h, t_step*s6, L, N, first_all_x);
}

void st8(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N, int first_all_x){
	const double s6 = 0.3595846493499923;

	st6(x, y, J, h, t_step*s6, L, N, first_all_x);
	st6(x, y, J, h, t_step*s6, L, N, first_all_x);
	st6(x, y, J, h, t_step*(1-4*s6), L, N, first_all_x);
	st6(x, y, J, h, t_step*s6, L, N, first_all_x);
	st6(x, y, J, h, t_step*s6, L, N, first_all_x);
}

void taylor17(double complex *x, double complex *y, double *J, double *h, double complex t_step, unsigned L, unsigned long N){
	double complex *z = y + N;

	memcpy(y, x, N * sizeof(double complex));

	for(unsigned i = 1; i <= 17; i++){
		const double complex norm = conj(t_step) / i;

		apply_hamiltonian(y, z, J, h, L, N);

#pragma omp parallel for
		for(unsigned long k = 0; k < N; k++){
			z[k] *= norm;
			x[k] += z[k];
		}

		double complex *dummy = y;
		y = z;
		z = dummy;
	}
}
