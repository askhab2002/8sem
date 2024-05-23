#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Fourier.hpp"
#define eps1 1e-5
#define eps2 1e-5

using namespace std;

double *Fourier_2d(function_pointer f, double p, int N) {

	double h = 1/((double)N - 0.5);

	double *F = (double *) calloc(N, sizeof(double));

        for(int i = 0; i < N; i++) {
                F[i] = (*f)(i * h);
        }

	cout << "   Функция или массив?(0, 1) " << endl;
	int index = 0;
	cin >> index;
	if(index == 1) {
		for(int i = 0; i < N/2; i++) {
			F[i] = i;
		}

		index = 0;
		for(int i = N - 1; i >= N/2; i--) {
			F[i] = -index;
			index++;
		}
	}

	double *D = coeff_out(F, N);

        double *C = (double *)calloc(N, sizeof(double));

        index = 0;
	for(int i = 0; i < N; i++) {
                if(i == 0 && fabs(Lambda(i, p, N)) < eps2) {
                       C[i] = 0;
		       index = 1;
                       continue;
                }
                C[i] = D[i]/Lambda(i, p, N);
        }

	double *Y = (double *)calloc(N, sizeof(double));

        for(int k = 0; k < N; k++) {
              Y[k] = fourier(C, N, k * h);

        }

	double *errors = error(Y, F, p, N);

	double error_norm = 0;
	double b_norm = 0;

        for(int i = 0; i < N; i++) {
                  error_norm += errors[i] * errors[i];
		  b_norm += F[i] * F[i];

        }

	double L_0 = 0;
	double L_0_ = 0;

	for(int i = 0; i < N; i++) {
		if(fabs(errors[i]) > L_0) {
			L_0 = fabs(errors[i]);
		}

		if(fabs(errors[i])/Y[i] > L_0_) {
			L_0_ = fabs(errors[i])/Y[i];
		}
	}

//	cout << "error_norm = " << sqrt(error_norm) << "   relative error = " << sqrt(error_norm)/sqrt(b_norm) << endl;

//	cout << " Fourier L_0 error = " << L_0 << "   Fourier L_0 relative = " << L_0_ << endl;
//        free(Y);
        free(D);
        free(C);
	free(errors);

	free(F);

	
	return Y;
}

double *error(double *Y, double *f, double p, int N) {
	double *errors = (double *)calloc(N, sizeof(double));

	errors[0] = f[0] - p * Y[0] - 2 * (N - 0.5) * (N - 0.5) * (Y[0] - Y[1]);
	errors[N - 1] = f[N - 1] - p * Y[N - 1] - (N - 0.5) * (N - 0.5) * (Y[N - 1] - Y[N - 2]);

	for(int i = 1; i < N - 1; i++) {
		errors[i] = f[i] - p * Y[i] + (N - 0.5) * (N - 0.5) * (Y[i + 1] - 2 * Y[i] + Y[i - 1]);
	}

	return errors;
}

double fourier(double *coeff, int num, double x) {

        double sum = 0;

        for(int i = 0; i < num ; i++) {
                sum += coeff[i] * cos(M_PI * i * x);
        }

        return sum;
}

double Lambda(int n, double p, int N) {
        return p + 2 * ((double)N - 0.5) * ((double)N - 0.5) * (1 - cos((M_PI * n)/((double)N - 0.5)));
}

double Scalar_(double *func, double h, int m, int n) {
        double scalar = 0;

	scalar += func[0] * cos(0) / 2;

	for(int k = 1; k < n; k++) {
		scalar += func[k] * cos(M_PI * k * h * m);
	}

	return scalar * h;
}

double *coeff_out(double *func, int values_number) {

	double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
	double h = 1/((double)values_number - 0.5);

	for(int m = 0; m < values_number ; m++) {
		scalar = Scalar_(func, h, m, values_number);

		coeff[m] = scalar * 2;

	}

	double sum = 0;
        for(int i = 1; i < values_number; i++) {
	        sum += coeff[i];
	}
        coeff[0] = func[0] - sum;



	return coeff;
}

