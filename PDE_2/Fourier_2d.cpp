#define _USE_MATH_DEFINES
#define eps 1e-6

#include <iostream>
#include <time.h>
#include <random>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "Fourier_2d.hpp"

using namespace std;

typedef double (* function_pointer)(double x);
typedef double (* function_pointer2)(double *coeff, int num, double x);
typedef double (* function_pointer_2d)(double x, double y);
typedef double (* function_pointer2_2d)(double **coeff, int num, double x, double y);


double **do_matrix(function_pointer_2d func, double *nodes, int node_number);
double fourier(double **coeff, int num, double x, double y);
double Scalar(double  **func, double h, int m, int n, int second);
double *generate_equel(double left, int num);
void coeff_out(double **func, int values_number, int second, double *coeff);
void coeff_out_2d(int node_number, double **func, double **coeff_new);
double Lambda_out(int n, double p, int N);



void coeff_out_2d(int node_number, double **func_matrix, double **coeff_new) {

        double **coeff = (double **) malloc(node_number * sizeof(double *));

        for(int i = 0; i < node_number; i++) {
		coeff[i] = (double *)calloc(node_number, sizeof(double));
                coeff_out(func_matrix, node_number, i, coeff[i]);
        }



        for(int i = 0; i < node_number; i++) {
                coeff_out(coeff, node_number, i, coeff_new[i]);
        }

	for(int i = 0; i < node_number; i++) {
                free(coeff[i]);
        }
        free(coeff);
        

        return;
}

void coeff_out(double **func, int values_number, int second, double *coeff) {

//        double *coeff = (double *) calloc(values_number , sizeof(double));
        double scalar = 0;
        double h = 1/((double)values_number - 0.5);


        for(int m = 0; m < values_number ; m++) {
                scalar = Scalar(func, h, m, values_number, second);

                coeff[m] = scalar * 2;
        }


        double sum = 0;
        for(int i = 1; i < values_number; i++) {
                sum += coeff[i];
        }
        coeff[0] = func[0][ second] - sum;

        return; //coeff;
}

double *generate_equel(double left, int num) {
        double step = 1/(num - 0.5);
        double *nodes = (double *)calloc( num, sizeof(double) );

        if(num >= 2) {
                nodes[0] = left;
        }

        for(int i = 1; i < num; i++) {
               nodes[i] = nodes[i - 1] + step;
        }

        return nodes;
}

double Scalar(double  **func, double h, int m, int n, int second) {

        double scalar = 0;

        scalar += func[0][second] * cos(0) / 2;

        for(int k = 1; k < n; k++) {
                scalar += func[k][second] * cos(M_PI * k * h * m);
        }


        return scalar * h;
}

double fourier(double **coeff, int num, double x, double y) {
	double sum = 0;

        for(int i = 0; i < num ; i++) {
		for(int j = 0; j < num; j++) {
			sum += coeff[i][j] * cos(M_PI * i * x) * cos(M_PI * j * y);
		}
        }

        return sum;
}

double **do_matrix(function_pointer_2d func, double *nodes, int node_number) {
	double **matrix = (double **) malloc( node_number * sizeof(double *));

	for(int i = 0; i < node_number; i++) {
		matrix[i] = (double *) malloc(node_number * sizeof(double));
		for(int j = 0; j < node_number; j++) {
			matrix[i][j] = (*func)(nodes[i], nodes[j]);
		}
	}

	return matrix;
}

double Lambda_out(int n, double p, int N) {
        return p + 2 * ((double)N - 0.5) * ((double)N - 0.5) * (1 - cos((M_PI * n)/((double)N - 0.5)));
}
