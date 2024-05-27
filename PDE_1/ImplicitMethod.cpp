#include <iostream>
#include <stdio.h>
#include <math.h>
#include "ImplicitMethod.hpp"

using namespace std;


void Sweep(double *y, double *a, double *b, double *c, double *f, int N) {
/*
        for(int i = 0; i < N - 1; i++) {
		a[i] = -a[i];
		b[i] = -b[i];
	}
*/
	double *alpha = (double *)calloc(N, sizeof(double));
	double *beta = (double *)calloc(N, sizeof(double));

	alpha[0] = b[0]/c[0];
	beta[0] = f[0]/c[0];

//	cout << " alpha_0 = " << alpha[0] << endl;
//	cout << " beta_0 = " << beta[0] << endl;

	for(int i = 1; i < N; i++) {
		alpha[i] = b[i]/(c[i] - a[i - 1] * alpha[i - 1]); 
		beta[i] = (f[i] + a[i - 1] * beta[i - 1])/(c[i] - a[i - 1] * alpha[i - 1]);
	}

//	cout << " alpha_N = " << alpha[N - 1] << endl;
//	cout << " beta_N = " << beta[N - 1] << endl;

//	double *y = (double *)calloc(N + 1, sizeof(double));

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
//	cout << "   Sweep norm = " << norm << endl;


	free(y_);


        return;
}


double **SweepMethod(function_pointer_2 f, function_pointer_1 p, function_pointer_1 u_0, int T, int N) {

	double *a = (double *)calloc(N - 1, sizeof(double));
	double *b = (double *)calloc(N - 1, sizeof(double));
	double *c = (double *)calloc(N, sizeof(double));

	double tau = 1/(double)T;
        double h = 1/((double)N - 0.5);

	double **u = (double **)malloc(T * sizeof(double *));
	for(int i = 0; i < T; i++) {
	       u[i] = (double *)calloc(N, sizeof(double));
        }	


	for(int i = 1; i < N - 2; i++) {
		a[i] = 1/(h * h);
		b[i] = 1/(h * h);
	}

	a[0] = 2/(h * h);
	a[N - 2] = 1/(h * h);
	b[0] = 1/(h * h);
	b[N - 2] = 1;

	for(int i = 0; i < N; i++) {
		u[0][i] = (*u_0)(i * h);
	}

	double *F = (double *)calloc(N, sizeof(double));

	for(int i = 1; i < T; i++) {
		for(int j = 0; j < N - 1; j++) {
			c[j] = (1/tau) + (2/(h * h)) + (*p)(j * h);
//			cout << c[j] << " ";
		}
		c[N - 1] = 1;
//		cout << c[N - 1] << endl;

		for(int j = 0; j < N - 1; j++) {
			F[j] = (u[i - 1][j]/tau) + (*f)(i * tau, j * h);
//			cout << F[j] << " ";
		}
//		cout << endl;
                
		//cout << " Tau = " << tau << " h = " << h << endl;
		Sweep(u[i], b, a, c, F, N - 1);
		u[i][N - 1] = u[i][N - 2];

		//cout << "a = " << a[0] << " b = " << b[0] << " c = " << c[0] << endl;;
	}

	free(a);
	free(b);
	free(c);
	free(F);

	return u;
}


			






