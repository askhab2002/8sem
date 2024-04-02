#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Fourier.hpp"
#include "SweepMethod.hpp"
#define eps1 1e-5

using namespace std;

double function1(double x) {

	return cos(M_PI * x );
}

double function2(double x) {
        if(fabs(x - 1) < eps1) {
                return 0;
        }

        if(fabs(x) < eps1) {
                return 0;
        }

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1));
}

double function3(double x) {

	double k = x;
	k = 1;

	return k;
}

double function4(double x) {

	if(x > 0.5) {
		return x - 1;
	}

	return x;
}


int main(void) {

	function_pointer f = function2;
	function_pointer p = function3;


	int N = 10;
	cout << "   Введите число узлов " << endl;
	cin >> N;

	double p_const = 1;

	double *x = (double *)calloc(N, sizeof(double));
	double *y_fourier = Fourier_2d(f, p_const, N);
	double *y_sweep = SweepMethod(f, p, N);

	FILE *file;
	FILE *file_1;
	FILE *file_2;

	double h = 1/((double)N - 0.5);

	file = fopen("RealY.txt", "w+");
	file_1 = fopen("FourierOut.txt", "w+");
	file_2 = fopen("SweepOut.txt", "w+");

	for(int i = 0; i < N; i++) {
                x[i] = i * h;
		fprintf(file_1, "%lf %lf\n", x[i], y_fourier[i]);
		fprintf(file_2, "%lf %lf\n", x[i], y_sweep[i]);
//		cout << i + 1 << "  Fourier: " << y_fourier[i] << "   Sweep: " << y_sweep[i] << endl;
	}

	double *x_ = (double *)calloc(5 * N, sizeof(double));
	double *y = (double *)calloc(5 * N, sizeof(double));
	for(int i = 0; i < 5 * N; i++) {
		x_[i] = i * (h/5);
		if(x_[i] > 1) {
			x_[i] = x_[i - 1];
		}

		y[i] = (*f)(x_[i]);
		fprintf(file, "%lf %lf\n", x_[i], y[i]);
	}

        free(x_);
	free(y);
        free(x);
	free(y_fourier);
	free(y_sweep);

	fclose(file);
	fclose(file_1);
	fclose(file_2);

	return 0;
}
