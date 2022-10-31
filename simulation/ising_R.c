#include <R.h>
#include <Rinternals.h>
#include <complex.h>

#include "ising_noisy_trace.h"

static R_INLINE double complex toC99(const Rcomplex *x){
#if __GNUC__
	double complex ans = (double complex) 0; /* -Wall */
	__real__ ans = x->r;
	__imag__ ans = x->i;
	return ans;
#else
	return x->r + x->i * I;
#endif
}

static R_INLINE void SET_C99_COMPLEX(Rcomplex *x, R_xlen_t i, double complex value){
	Rcomplex *ans = x+i;
	ans->r = creal(value);
	ans->i = cimag(value);
}

SEXP noisy_trace(SEXP L, SEXP J, SEXP h, SEXP t, SEXP t_step, SEXP n_sources, SEXP scheme, SEXP time_series, SEXP full_data, SEXP operator_id, SEXP first_all_x, SEXP seed){
	const unsigned long n_steps = (unsigned long) ceil(cabs(toC99(COMPLEX(t))) / asReal(t_step));
	const unsigned long n = (asInteger(time_series)? n_steps:1L) * (asInteger(full_data)? asInteger(n_sources):2);
	SEXP out = PROTECT(allocVector(CPLXSXP, n)); // need space for rho

	double *coupling = malloc(3*sizeof(double));
	memcpy(coupling, REAL(J), 3*sizeof(double));
	double *field = malloc(asInteger(L)*sizeof(double));
	memcpy(field, REAL(h), asInteger(L)*sizeof(double));

	double complex *trace;
	trace = trace_estimator(coupling, field, toC99(COMPLEX(t)), asInteger(L), n_steps, asInteger(n_sources), asInteger(scheme), asInteger(time_series), asInteger(full_data), asInteger(operator_id), asInteger(first_all_x), asInteger(seed));

	for(unsigned long i = 0; i < n; i++) SET_C99_COMPLEX(COMPLEX(out), i, trace[i]);

	free(coupling);
	free(field);
	free(trace);
	UNPROTECT(1);

	return out;
}
