#include <stdio.h>
#include <iostream>
#include "Coefficients.hpp"
#include "Functions.hpp"

using namespace std;

typedef void  (* ans_func)(double x, double *y, int n);

void test(double x_0, double y_0, int steps, int n, function_pointer f, ans_func g, FILE *file);

int main(void) {

	int n = 1;
	int steps = 10;
	double x_0 = 0;
	double y_0 = 0;
        
	FILE *file1 = fopen("res1.txt", "w+");
	FILE *file2 = fopen("res2.txt", "w+");
	FILE *file3 = fopen("res3.txt", "w+");
	FILE *file4 = fopen("res4.txt", "w+");
	FILE *file6 = fopen("res6.txt", "w+");

	function_pointer F = F1;
        ans_func G = G1;

	test(x_0, y_0, steps, n, F, G, file1);

	F = F2;
        G = G2;

        test(x_0, y_0, steps, n, F, G, file2);

	F = F3;
        G = G3;

        test(x_0, y_0, steps, n, F, G, file3);

	F = F4;
	G = G4;

	test(x_0, y_0, steps, n, F, G, file4);

	F = F6;
	G = G6;
	y_0 = 1;
	test(x_0, y_0, steps, n, F, G, file6);


	return 0;

}	

void test(double x_0, double y_0, int steps, int n, function_pointer F, ans_func g, FILE *file) {
	double h = 1/((double)steps);
        double *y = (double *)calloc(n, sizeof(double));
        double *y_n = (double *)calloc(n, sizeof(double));
        double *E_n = (double *)calloc(n, sizeof(double));

        double *k1 = (double *)calloc(n, sizeof(double));
        double *k2 = (double *)calloc(n, sizeof(double));
        double *k3 = (double *)calloc(n, sizeof(double));
        double *k4 = (double *)calloc(n, sizeof(double));

	double *real_y = (double *)calloc(n, sizeof(double));

	for(int i = 0; i < n; i++) {
		y[i] = y_0;
	}
	double x = x_0;

	for(int i = 0; i < steps; i++) {
                Yn(F, x, y, h, n, y_n, k1, k2, k3, k4);
                En(F, x, y, h, n, E_n, k1, k2, k3, k4);

                x = x + h;
                for(int j = 0; j < n; j++) {
                        y[j] = y_n[j];

                }


                fprintf(file, "x = %lf ", x);

		g(x, real_y, n);

		fprintf(file, "real_y = ");

		for(int j = 0; j < n; j++) {
			fprintf(file, "%lf ", real_y[j]);
		}

		fprintf(file, "approx_y = ");

		for(int j = 0; j < n; j++) {
                        fprintf(file, "%lf ", y[j]);
                }

		

		fprintf(file, "Norm|y() - yh| = %lf ", L0(y_n, real_y, n));
		fprintf(file, "Norm|E| = %lf ", L0(E_n, n));

		fprintf(file, "\n\n");

                cout << y_n[0] << " ";

        }

	cout << endl;


        fclose(file);
        free(y);
	free(real_y);
        free(y_n);
	free(E_n);
        free(k1);
        free(k2);
        free(k3);
        free(k4);

	return;
}

