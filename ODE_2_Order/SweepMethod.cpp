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

	cout << " alpha_0 = " << alpha[0] << endl;
	cout << " beta_0 = " << beta[0] << endl;

	for(int i = 1; i < N; i++) {
		alpha[i] = b[i]/(c[i] - a[i - 1] * alpha[i - 1]); 
		beta[i] = (f[i] + a[i - 1] * beta[i - 1])/(c[i] - a[i - 1] * alpha[i - 1]);
	}

	cout << " alpha_N = " << alpha[N - 1] << endl;
	cout << " beta_N = " << beta[N - 1] << endl;

	double *y = (double *)calloc(N + 1, sizeof(double));

	y[N] = (f[N] + a[N - 1] * beta[N - 1])/(c[N] - a[N - 1] * alpha[N - 1]);

	for(int i = N - 1; i >= 0; i--) {
	       y[i] = alpha[i] * y[i + 1] + beta[i];
        }

	free(alpha);
	free(beta);

	double *y_ = (double *)calloc(N + 1, sizeof(double));

	y_[0] = c[0] * y[0] + b[0] * y[1];
	y_[N] = a[N - 1] * y[N - 1] + c[N] * y[N];
	for(int i = 1; i < N; i++) {
		y_[i] = a[i - 1] * y[i - 1] + c[i] * y[i] + b[i] * y[i + 1];
	}

	double norm = 0;
	for(int i = 0; i < N + 1; i++) {
		norm += (y_[i] - f[i]) * (y_[i] - f[i]);
	}

	norm = sqrt(norm/((double)N + 0.5));
	cout << "   Sweep norm = " << norm << endl;

	free(y_);


        return y;
}

double *SweepMethod(function_pointer f, function_pointer p, int N) {
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
		a[i] = 1/(h * h);
		b[i] = 1/(h * h);
	}

	b[0] = 2/(h * h);
	c[0] = 2/(h * h) + P[0];
	a[N - 2] = 1/(h * h);
	c[N - 1] = 1/(h * h) + P[N - 1];
//	F[N - 1] = 0;

	double *y =  Sweep(a, b, c, F, N - 1);

	free(a);
	free(b);
	free(c);
	free(F);
	free(P);
	

	return y;
}






