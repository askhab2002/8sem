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

	double *y_fourier = Fourier_2d(f, p_const, N);
	double *y_sweep = SweepMethod(p, f, N);

	for(int i = 0; i < N; i++) {
//		cout << i + 1 << "  Fourier: " << y_fourier[i] << "   Sweep: " << y_sweep[i] << endl;
	}

	free(y_fourier);
	free(y_sweep);

	return 0;
}
