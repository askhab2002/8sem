#include <iostream>
#include <stdio.h>
#include <math.h>
#include "SweepMethod.hpp"

using namespace std;

double *Sweep(double *a, double *b, double *c, double *f, int N) {
	double *alpha = (double *)calloc(N, sizeof(double));
	double *beta = (double *)calloc(N, sizeof(double));

	alpha[0] = b[0]/c[0];
	beta[0] = f[0]/c[0];

	for(int i = 1; i < N; i++) {
		alpha[i] = b[i]/(c[i] - a[i - 1] * alpha[i - 1]); 
		beta[i] = (f[i] + a[i - 1] * alpha[i - 1])/(c[i] - a[i - 1] * alpha[i - 1]);
	}

	double *y = (double *)calloc(N + 1, sizeof(double));

	y[N] = (f[N] + a[N - 1] * beta[N - 1])/(c[N] - a[N - 1] * alpha[N - 1]);

	for(int i = N - 1; i >= 0; i--) {
	       y[i] = alpha[i] * y[i + 1] + beta[i];
        }

	free(alpha);
	free(beta);

        return y;
}

double *SweepMethod(function_pointer p, function_pointer f, int N) {
	double *F = (double *)calloc(N, sizeof(double));
	double *P = (double *)calloc(N, sizeof(double));

	double h = 1/((double)N - 0.5);

	for(int i = 0; i < N; i++) {
		F[i] = (*f)(i * h);
		P[i] = (*p)(i * h);
	}

	double *a = (double *)calloc(N - 1, sizeof(double));
	double *b = (double *)calloc(N - 1, sizeof(double));
	double *c = (double *)calloc(N, sizeof(double));

	for(int i = 0; i < N; i++) {
		c[i] = 2/(h * h) + P[i];
	}

	for(int i = 0; i < N - 1; i++) {
		a[i] = -1/(h * h);
		b[i] = -1/(h * h);
	}
	

	double *y =  Sweep(a, b, c, F, N - 1);

	free(a);
	free(b);
	free(c);
	free(F);
	free(P);
	

	return y;
}






