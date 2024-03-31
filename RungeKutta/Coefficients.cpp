#include <math.h>
#include <iostream>

using namespace std;


typedef void  (* function_pointer)(double x, double *y, int n, double *temporary);


void K(function_pointer f, double x, double *y, double h, int n, double *k_1) {
	f(x, y, n, k_1);

	for(int i = 0; i < n; i++) {
		k_1[i] = h * k_1[i];
	}

	return;
}



void Yn(function_pointer f, double x, double *y, double h, int n, double *y_n, double *k1, double *k2, double *k3, double *k4) {
	K(f, x, y, h, n, k1);

        for(int i = 0; i < n; i++) {
		y_n[i] = y[i] + 0.5 * k1[i];
	}

	K(f, x + 0.5 * h, y_n, h, n, k2);

	for(int i = 0; i < n; i++) {
                y_n[i] = y[i] + 0.5 * k2[i];
        }

	K(f, x + 0.5 * h, y_n, h, n, k3);

	for(int i = 0; i < n; i++) {
                y_n[i] = y[i] + k3[i];
        }

        K(f, x + h, y_n, h, n, k4);

	for(int i = 0; i < n; i++) {
		y_n[i] = y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])/6;
	}

//	cout << " 0 = " << y_n[0] << " ";

	return;
}

void En(function_pointer f, double x, double *y, double h, int n, double *E_n, double *k1, double *k2, double *k3, double *k4) {
        K(f, x, y, h, n, k1);

        for(int i = 0; i < n; i++) {
                E_n[i] = y[i] + 0.5 * k1[i];
        }

        K(f, x + 0.5 * h, E_n, h, n, k2);

        for(int i = 0; i < n; i++) {
                E_n[i] = y[i] + 0.5 * k2[i];
        }

        K(f, x + 0.5 * h, E_n, h, n, k3);

        for(int i = 0; i < n; i++) {
                E_n[i] = y[i] + k3[i];
        }

        K(f, x + h, E_n, h, n, k4);

        for(int i = 0; i < n; i++) {
                E_n[i] = (((double)2)/(3)) * (k1[i] - k2[i] - k3[i] + k4[i]);
	}

        return;
}


double L2_h(double *x, double *y, double h, int n) {
	double norm = 0;

	for(int i = 0; i < n; i++) {
		norm = norm + x[i] * y[i];
	}

	return h * norm;
}

double L0(double *x, int n) {

	double max = 0;

	for(int i = 0; i < n; i++) {
		if(max < fabs(x[i])) {
			max = fabs(x[i]);
		}
	}

	return max;
}

double L0(double *x, double *y, int n) {

        double max = 0;

        for(int i = 0; i < n; i++) {
                if(max < fabs(x[i] - y[i])) {
                        max = fabs(x[i] - y[i]);
                }
        }

        return max;
}

double L2_h_diff(double *x, double *y, double h, int n) {
	double norm = 0;

        for(int i = 0; i < n; i++) {
                norm = norm + (x[i] - y[i]) * (x[i] - y[i]);
        }

        return h * norm;
}
