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

	return cos(M_PI * x)/(M_PI * M_PI + 1);
}

double function5(double x) {

        return cos(x);
}

double function6(double x) {

        return 0.5 * (cos(x) + (exp(1 - x) * (exp(2 * x) + 1) * sin(1))/(exp(2) - 1));
}

double function7(double x) {

	return exp(x);
}

double function8(double x) {

	return 0.5 * exp(x) * (x - 1);
}

function_pointer choose_func(int k) {

	switch(k) {
		case 1: 
			return function1;
		case 2:
                        return function2;
		case 3:
                        return function3;
		case 4:
                        return function4;
		case 5:
                        return function5;
                case 6:
                        return function6;
                case 7:
                        return function7;
                case 8:
                        return function8;
	}

	return function1;

}



void test(int sol_n, int f_n, int p_n) {

	function_pointer solution = choose_func(sol_n);
	function_pointer f = choose_func(f_n);
	function_pointer p = choose_func(p_n);


	int N = 10;
	cout << "   Введите число узлов " << endl;
	cin >> N;

	double p_const = (*p)(0.5);

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

	double *x_ = (double *)calloc(N, sizeof(double));
	double *y = (double *)calloc(N, sizeof(double));

	double L_0_error_Fourier = 0;
	double L_0_error_Sweep = 0;

	for(int i = 0; i < N; i++) {
		x_[i] = i * (h);
		if(x_[i] > 1) {
			x_[i] = x_[i - 1];
		}

		y[i] = (*solution)(x_[i]);

		if(L_0_error_Fourier < fabs(y[i] - y_fourier[i])) {
			L_0_error_Fourier = fabs(y[i] - y_fourier[i]);
		}

		if(L_0_error_Sweep < fabs(y[i] - y_sweep[i])) {
                        L_0_error_Sweep = fabs(y[i] - y_sweep[i]);
                }
		fprintf(file, "%lf %lf\n", x_[i], y[i]);
	}

	cout << "L_0 error Fourier = " << L_0_error_Fourier << endl;
        cout << "L_0 error Sweep = " << L_0_error_Sweep << endl;


        free(x_);
	free(y);
        free(x);
	free(y_fourier);
	free(y_sweep);

	fclose(file);
	fclose(file_1);
	fclose(file_2);
}

void test_sweep(int sol_n, int f_n, int p_n) {

        function_pointer solution = choose_func(sol_n);
        function_pointer f = choose_func(f_n);
        function_pointer p = choose_func(p_n);


        int N = 10;
        cout << "   Введите число узлов " << endl;
        cin >> N;

        double *x = (double *)calloc(N, sizeof(double));
        double *y_sweep = SweepMethod(f, p, N);

        FILE *file;
        FILE *file_2;

        double h = 1/((double)N - 0.5);

        file = fopen("RealY.txt", "w+");
        file_2 = fopen("SweepOut.txt", "w+");

	for(int i = 0; i < N; i++) {
                x[i] = i * h;
		fprintf(file_2, "%lf %lf\n", x[i], y_sweep[i]);
//		cout << i + 1 << "  Fourier: " << y_fourier[i] << "   Sweep: " << y_sweep[i] << endl;
	}

	double *x_ = (double *)calloc(N, sizeof(double));
	double *y = (double *)calloc(N, sizeof(double));

	double L_0_error_Sweep = 0;

	for(int i = 0; i < N; i++) {
		x_[i] = i * (h);
		if(x_[i] > 1) {
			x_[i] = x_[i - 1];
		}

		y[i] = (*solution)(x_[i]);


		if(L_0_error_Sweep < fabs(y[i] - y_sweep[i])) {
                        L_0_error_Sweep = fabs(y[i] - y_sweep[i]);
                }
		fprintf(file, "%lf %lf\n", x_[i], y[i]);
	}

        cout << "L_0 error Sweep = " << L_0_error_Sweep << endl;


        free(x_);
	free(y);
        free(x);
	free(y_sweep);

	fclose(file);
	fclose(file_2);
}

int main(void) {

	test(6, 5, 3);
	//test_sweep(4, 1, 3);

	return 0;

}
