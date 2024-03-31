#include <math.h>

void F1(double x, double *y, int n, double *F) {
	double a  = x;
        x = y[0];
        x = a;

	for(int i = 0; i < n; i++) {
		F[i] = x;
	}
}

void G1(double x, double *y, int n) {
	for(int i = 0; i < n; i++) {
		y[i] = 0.5 * x * x;
	}
}

void F2(double x, double *y, int n, double *F) {
	double a  = x;
	x = y[0];
	x = a;

	for(int i = 0; i < n; i++) {
		F[i] = x * x;
	}
}

void G2(double x, double *y, int n) {
        for(int i = 0; i < n; i++) {
                y[i] = (x * x * x)/3;
        }
}


void F3(double x, double *y, int n, double *F) {
	double a  = x;
        x = y[0];
        x = a;

	for(int i = 0; i < n; i++) {
                F[i] = x * x * x;
        }
}

void G3(double x, double *y, int n) {
        for(int i = 0; i < n; i++) {
                y[i] = 0.25 * x * x * x * x;
        }
}


void F4(double x, double *y, int n, double *F) {
	double a  = x;
        x = y[0];
        x = a;

	for(int i = 0; i < n; i++) {
                F[i] = x * x * x * x;
        }
}

void G4(double x, double *y, int n) {
        for(int i = 0; i < n; i++) {
                y[i] = 0.2 * x * x * x * x * x;
        }
}


void F5(double x, double *y, int n, double *F) {
        double a  = x;
        x = y[0];
        x = a;

        for(int i = 0; i < n; i++) {
                F[i] = exp(x);
        }
}

void F6(double x, double *y, int n, double *F) {
        double a  = x;
        x = y[0];
        x = a;

        for(int i = 0; i < n; i++) {
                F[i] = y[i];
        }
}

void G6(double x, double *y, int n) {
        for(int i = 0; i < n; i++) {
                y[i] = exp(x);
        }
}


