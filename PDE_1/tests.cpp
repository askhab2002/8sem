#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
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

double function8(double x){

	return cos(M_PI * x) + cos(M_PI * 2 * x) + cos(M_PI * 3 * x) + cos(M_PI * 4 * x);
}

double function9(double t, double x) {

	double sum1 = cos(M_PI * x) * exp(-M_PI * M_PI * t) + cos(M_PI * 2 * x) * exp(-M_PI * M_PI * 4 * t);
	double sum2 = cos(M_PI * 3 * x) * exp(-M_PI * M_PI * 9 * t) + cos(M_PI * 4 * x) * exp(-M_PI * M_PI * 16 * t);

	return sum1 + sum2;
}

void result_out(double **u, double **u_real, int N, int T, string filename) {

	fstream out;
        out.open(filename, std::ofstream::out | std::ofstream::trunc);

	double *x_grid = (double *)calloc(N, sizeof(double));
	double *t_grid = (double *)calloc(T, sizeof(double));

	double h = 1 / ((double)N - 0.5);
	double t = 1 / ((double)T);

	for(int i = 0; i < N; i++) {
		x_grid[i] = i * h;
	}

	for(int i = 0; i < T; i++) {
		t_grid[i] = i * t;
	}

	for(int i = 0; i < N; i++) {

		for(int j = 0; j < T; j++) {
			out << x_grid[i] << " " << t_grid[j] << " " << u[j][i] << " " << u_real[j][i] << endl;

		}
	}

	free(x_grid);
	free(t_grid);

	out.close();
}

function_pointer_1 choose_func(int k) {

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
                case 8:
                        return function8;
	}

	return function1;

}

function_pointer_2 choose_func_2(int k) {

        switch(k) {
                case 6:
                        return function6;
		case 7:
			return function7;
                case 9:
                        return function9;
        }

        return function6;

}



void test(int func1, int func2, int func3, int func4, int T, int N) {

	function_pointer_1 u_0 = choose_func(func1); //function5;
	function_pointer_1 p = choose_func(func2); //function4;
	function_pointer_2 f = choose_func_2(func3); //function7;


	double **u = SweepMethod(f, p, u_0, T, N);

	function_pointer_2 u_real_ = choose_func_2(func4); //function6

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

			error_mean += fabs(u[i][j] - u_real[i][j]);

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

	double mean_error = error_mean/(N * T);

	cout << " L_0 error = " << L_0_error << " L_0 relative error = " << L_0_relative_error << endl;
        cout << " error mean disp = " << mean_error << endl;

	
	result_out(u, u_real, N, T, "result.txt");

	for(int i = 0; i < T; i++) {
		free(u[i]);
		free(u_real[i]);
	}
	free(u);
	free(u_real);
	

	return;
}

int main(void) {
	int T = 10;
	int N = 10;

	cout << " Введите число разбиения пространственной координаты" << endl;
        cin >> N;

        cout << " Введите число разбиения временной координаты" << endl;
        cin >> T;

	int func_u_0 = 8; //5;
	int func_p = 4;
	int func_f = 7; //7;
	int func_solution = 9; //6;

	//test(func_u_0, func_p, func_f, func_solution, T, N);
	
	func_u_0 = 5;
        func_p = 4;
        func_f = 7;
        func_solution = 6;

        test(func_u_0, func_p, func_f, func_solution, T, N);
	

	return 0;
}
