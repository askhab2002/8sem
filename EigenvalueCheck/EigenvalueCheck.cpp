#include <iostream>
#include <cmath>
#include <math.h>
#include "EigenvalueCheck.hpp"

using namespace std;

double Y_k(int k, int n, int N) {
       return cos((M_PI	* n * k) / ((double)N - 0.5));
}

double lambda_n(double p, int n, int N) {
	return - p + 2 * (N - 0.5) * (N - 0.5) * (1 - cos((M_PI * n)/(N - 0.5)));
}


int main(void) {
	cout << "  Введите число N" << endl;

	double N = 0;
	cin >> N;

	cout << "  Введите число p" << endl;

	double p = 0;
	cin >> p;

	double h = 1/((double)N - 0.5);

	double lambda = 0;
	double value_right = 0;
	double value_left = 0;

	for(int n = 0; n < N; n++) {
		lambda = lambda_n(p, n, N);
   
		for(int k = 1; k < N; k++) {
			value_right = (Y_k(k + 1, n, N) - 2 * Y_k(k, n, N) + Y_k(k - 1, n, N))/(h * h) + p * Y_k(k, n, N);
			value_left = -lambda * Y_k(k, n, N);

//                        cout << value_right << " " << value_left << endl;
			if(fabs(value_right - value_left) > 1e-5) {
				
				cout << "   Не сходится, всё провалилось!" << endl;
				return 0;
				break;
			}
		}
	}

	cout << "  Собственное значение сошлось!" << endl;

	return 0;
}



