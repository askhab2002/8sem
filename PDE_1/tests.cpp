#include <iostream>
#include <stdio.h>
#include <math.h>
#include "ImplicitMethod.hpp"
#define eps 1e-8

using namespace std;

double function1(double x) {

	return cos(M_PI * x );
}

double function2(double x) {
        if(fabs(x - 1) < eps) {
                return 0;
        }

        if(fabs(x) < eps) {
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

        double k = x;
        k = 0;

        return k;
}

double function5(double x) {

        return cos(M_PI * x);
}


double function6(double t, double x) {

        return cos(M_PI * x) * exp(-M_PI * M_PI * t);
}

double function7(double x, double t) {
        double k = x;
	k = t;
	k = 0;

	return k;
}	



int main(void) {

	function_pointer_1 u_0 = function5;
	function_pointer_1 p = function4;
	function_pointer_2 f = function7;

	int T = 3;
        int N = 10;

	cout << " Введите число разбиения пространственной координаты" << endl;
	cin >> N;

	cout << " Введите число разбиения временной координаты" << endl;
	cin >> T;

	double **u = SweepMethod(f, p, u_0, T, N);

	function_pointer_2 u_real_ = function6;

	double h = 1/((double)N - 0.5);
	double t = 1/(double)T;

	double **u_real = (double **)malloc(T * sizeof(double *));
        for(int i = 0; i < T; i++) {
		u_real[i] = (double *)calloc(N, sizeof(double *));
		for(int j = 0; j < N; j++) {
			u_real[i][j] = (*u_real_)(i * t, j * h);
		}
	}

        int index = 0;
	cout << " Нужно ли выводить функцию на сетка?(1 - да, 0 - нет)" << endl;
	cin >> index;

	if(index == 1) {

	        cout << " Real: " << endl;
                for(int i = 0; i < T; i++) {
		        for(int j = 0; j < N; j++) {
			        cout << u_real[i][j] << " ";
		        }
		        cout << endl;
	        }


                cout << " Numerical: " << endl;
                for(int i = 0; i < T; i++) {
                        for(int j = 0; j < N; j++) {
                                cout << u[i][j] << " ";
                        }
                        cout << endl;
		}
	}

	double L_0_error = 0;
	double L_0_relative_error = 0;
	double error_mean = 0;

	for(int i = 0; i < T; i++) {
		for(int j = 0; j < N - 1; j++) {

			error_mean += fabs(u[i][j] - u_real[i][j]) * fabs(u[i][j] - u_real[i][j]);

			if(L_0_error < fabs(u[i][j] - u_real[i][j])) {
				L_0_error = fabs(u[i][j] - u_real[i][j]);
			}
                        
			if(fabs(u[i][j]) < eps) {
				continue;
			}
			if(L_0_relative_error < fabs((u[i][j] - u_real[i][j]) / u[i][j])) {
                                L_0_relative_error = fabs((u[i][j] - u_real[i][j]) / u[i][j]);
                        }
		}
	}

	cout << " L_0 error = " << L_0_error << " L_0 relative error = " << L_0_relative_error << endl;
        cout << " error mean disp = " << sqrt(error_mean) << endl;


	for(int i = 0; i < T; i++) {
		free(u[i]);
		free(u_real[i]);
	}
	free(u);
	free(u_real);
	

	return 0;
}
