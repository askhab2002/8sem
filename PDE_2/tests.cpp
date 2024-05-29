#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include "BSolver.hpp"
#include "Fourier_2d.hpp"
#define eps 1e-6
#define eps1 1e-6


using namespace std;

double function1(double x);
double function2(double x);

double function1_2d(double x, double y);
double function2_2d(double x, double y);

typedef double (* function_pointer_3d)(double x, double y, double t);

double function1_3d(double x, double y, double t);

double function3_2d(double x, double y);


void StepMatrix(double **f_, function_pointer_3d f_3d, double **u_prev, int N_x, int N_y, int T, int k);
void fourier_matrix(double **coeff, int N_x, int N_y, double **u);

void step_solution(double **u_coeff, double **f_coeff, double **u_0, double **u, double **f_, function_pointer_3d f_3d, int N_x, int N_y, int T, int k);
void Fourier_solution(double **u_coeff, double **f_coeff, double **u_0, double **u, double **f_, function_pointer_3d f_3d, int N_x, int N_y, int T);

void print_sol(function_pointer_3d sol, int N_x, int N_y, int T);
void sol(function_pointer_3d sol, int N_x, int N_y, int T, int k);

double solution_3(double x, double y, double t);
double solution_1(double x, double y, double t);

double **Inverse_(int N_x, int N_y, int T);
double **TriangularSystemP_(int N_x, int N_y, int T, double *f);

int main(void) {
/*
	cout << "    Введите число корней" << endl;

	int N = 0;
	int mIter = 10;
	cin >> N;

	function_pointer f = function1;

	cout << "   Метод с обуславливателем  " << endl;
	test_3(N, f, mIter);
*/



	cout << " Тест на 2д Фурье решение уравнения теплопроводности" << endl;
        cout << endl;

	function_pointer_3d f_3d = function1_3d;
	function_pointer_2d u_0_ = function3_2d;

	int N_x = 5;
	int N_y = N_x;
	int T = 3;

	cout << " Введите число пространственного разбиения" << endl;
	cin >> N_x;
	N_y = N_x;

	cout << " Введите число временного разбиения" << endl;
	cin >> T;

	double h_x = 1/((double)N_x - 0.5);
        double h_y = 1/((double)N_y - 0.5);

	double **u_0 = (double **)malloc(N_x * sizeof(double *));

	cout << " Значение функции на шаге 0" << endl;
	cout << endl;

	for(int i = 0; i < N_x; i++) {
		u_0[i] = (double *)calloc(N_y, sizeof(double));
		for(int j = 0; j < N_y; j++) {
			u_0[i][j] = (*u_0_)(i * h_x, j * h_y);
			cout << u_0[i][j] << " ";
		}
		cout << endl;
		cout << endl;
	}
	cout << endl;


	double **u_coeff = (double **)malloc(N_x * sizeof(double *));
	double **f_coeff = (double **)malloc(N_x * sizeof(double *));
	double **u = (double **)malloc(N_x * sizeof(double *));
	double **f_ = (double **)malloc(N_x * sizeof(double *));

        for(int i = 0; i < N_x; i++) {
                u_coeff[i] = (double *)calloc(N_y, sizeof(double));
		f_coeff[i] = (double *)calloc(N_y, sizeof(double));
		u[i] = (double *)calloc(N_y, sizeof(double));
		for(int j = 0; j < N_y; j++) {
			u[i][j] = u_0[i][j];
		}
		f_[i] = (double *)calloc(N_y, sizeof(double));
	}



	Fourier_solution(u_coeff, f_coeff, u_0, u, f_, f_3d, N_x, N_y, T);

	function_pointer_3d sol = solution_3;
	print_sol(sol, N_x, N_y, T);

	for(int i = 0; i < N_x; i++) {
		free(u_coeff[i]);
		free(f_coeff[i]);
		free(u_0[i]);
		free(u[i]);
		free(f_[i]);
	}

	free(u_coeff);
        free(f_coeff);
        free(u_0);
        free(u);
        free(f_);


	return 0;
}




void sol(function_pointer_3d sol, int N_x, int N_y, int T, int k) {

	double h_x = 1/((double)N_x - 0.5);
        double h_y = 1/((double)N_y - 0.5);
        double tau = 1/(double)T;
        
	for(int i = 0; i < N_x; i++) {
		for(int j = 0; j < N_y; j++) {
			cout << (*sol)(i * h_x, j * h_y, tau * k) << " ";
		}
		cout << endl;
	}

	return;
}

void print_sol(function_pointer_3d f, int N_x, int N_y, int T) {

	
	for(int k = 0; k < T; k++) {
	       sol(f, N_x, N_y, T, k);
	       cout << endl;
        }

        return;
}	

	

void StepMatrix(double **f_, function_pointer_3d f_3d, double **u_prev, int N_x, int N_y, int T, int k) {
	double tau = 1/(double)T;
        double h_x = 1/((double)N_x - 0.5);
        double h_y = 1/((double)N_y - 0.5);

        for(int i = 0; i < N_x; i++) {

                for(int j = 0; j < N_y; j++) {
                        f_[i][j] = (*f_3d)(i * h_x, j * h_y, k * tau) + u_prev[i][j] / tau;
		}
        }

	return;
}

void fourier_matrix(double **coeff, int N_x, int N_y, double **u) {

	double h_x = 1/((double)N_x - 0.5);
        double h_y = 1/((double)N_y - 0.5);

	for(int i = 0; i < N_x; i++) {
		for(int j = 0; j < N_y; j++) {
			u[i][j] = fourier(coeff, N_x, i * h_x, j * h_y);
		}
	}

	return;
}

void step_solution(double **u_coeff, double **f_coeff, double **u_0, double **u, double **f_, function_pointer_3d f_3d, int N_x, int N_y, int T, int k) {
 
        //cout << " u_0: " << endl;
	for(int i = 0; i < N_x; i++) {
                for(int j = 0; j < N_y; j++) {
                        
                        u_0[i][j] = u[i][j];
			//cout << u_0[i][j] << " ";
                }
                //cout << endl;
        }

	StepMatrix(f_, f_3d, u_0, N_x, N_y, T, k);

	
	double *lambda_x = (double *)calloc(N_x, sizeof(double));
	double *lambda_y = (double *)calloc(N_y, sizeof(double));

	for(int i = 0; i < N_x; i++) {
		lambda_x[i] = Lambda_out(i, 0, N_x);
	}

	for(int i = 0; i < N_y; i++) {
                lambda_y[i] = Lambda_out(i, 0, N_y);
        }

	cout << " f_: " << endl;

        for(int i = 0; i < N_x; i++) {
                for(int j = 0; j < N_y; j++) {
                       cout << f_[i][j] << " ";
                }
                cout << endl;
        }

	coeff_out_2d(N_x, f_, f_coeff);
/*
	cout << " u_0: " << endl;
	for(int i = 0; i < N_x; i++) {
		for(int j = 0; j < N_y; j++) {
			cout << u_0[i][j] << " ";
			u_0[i][j] = u[i][j];
		}
		cout << endl;
	} */

	cout << " f_coeff: " << endl;
	for(int i = 0; i < N_x; i++) {
		for(int j = 0; j < N_y; j++) {
			cout << f_coeff[i][j] << " ";
			u_coeff[i][j] = f_coeff[i][j] / (lambda_x[i] + lambda_y[j] + T);
		}
		cout << endl;
	}

        fourier_matrix(u_coeff, N_x, N_y, u);


	free(lambda_x);
	free(lambda_y);

	return;
}

void Fourier_solution(double **u_coeff, double **f_coeff, double **u_0, double **u, double **f_, function_pointer_3d f_3d, int N_x, int N_y, int T) {
/*
	cout << endl;
        for(int i = 0; i < N_x; i++) {
                for(int j = 0; j < N_y; j++) {
                        cout << u_0[i][j] << " ";
                }

                cout << endl;
                cout << endl;
	}
*/

	for(int t = 1; t < T; t++) {
		step_solution(u_coeff, f_coeff, u_0, u, f_, f_3d, N_x, N_y, T, t);

		cout << " Значение функции на шаге " << t << endl;
		cout << endl;
		for(int i = 0; i < N_x; i++) {
                        for(int j = 0; j < N_y; j++) {
                                cout << u[i][j] << " ";
                        }
                        cout << endl;
		        cout << endl;
                }

	}

	return;
}


	

double function1_3d(double x, double y, double t) {
	x = 0;
	y = x;
	t = x;

	return x + y + t;
}

double function3_2d(double x, double y) {

	return cos(1 * M_PI * x) * cos(1 * M_PI * y);
}

double solution_3(double x, double y, double t) {

	return exp(- M_PI * M_PI * (1 + 1) * t) * cos(1 * M_PI * x) * cos(1 * M_PI * y);
}


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

double function1_2d(double x , double y) {
        return cos(x * 4 * M_PI) *  cos(y * 2 * M_PI);
}

double solution_1(double x, double y, double t) {

        return exp(- M_PI * M_PI * (16 + 4) * t) * cos(4 * M_PI * x) * cos(2 * M_PI * y);
}


double function2_2d(double x, double y) {
        y = 0;
        x += y;

        if(fabs(x - 1) < eps) {
                return 0;
        }

        if(fabs(x) < eps) {
                return 0;
        }

        return exp(1/((2 * x - 1) * (2 * x - 1) - 1)) ;
}
