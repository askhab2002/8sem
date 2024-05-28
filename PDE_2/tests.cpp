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

int main(void) {

	cout << "    Введите число корней" << endl;

	int N = 0;
	int mIter = 10;
	cin >> N;

	function_pointer f = function1;

	cout << "   Метод с обуславливателем  " << endl;
	test_3(N, f, mIter);



	function_pointer_2d func = function1_2d;
	int node_number = 6;

	double **coeff_new = coeff_out_2d(node_number, func);

	cout << " Матрица коэффициентов: " << endl;

	for(int i = 0; i < node_number; i++) {
		for(int j = 0; j < node_number; j++) {
			cout << coeff_new[i][j] << " ";
		}
		cout << endl;
	}

	return 0;
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
