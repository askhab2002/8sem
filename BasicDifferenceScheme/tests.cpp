#include <stdio.h>
#include <iostream>
#include "methods.hpp"

using namespace std;

int main(void) {

	double A_ = 1;
	int n = 1;
	int number_method = 1;

	cout << "   Введите значение параметра А" << endl;
	cin >> A_;

	cout << "   Введите значения порядка дробления узлов n" << endl;
	cin >> n;

	cout << "   Введите номер разностной схемы" << endl;
	cin >> number_method;

	cout << E_n(number_method, n, A_) << endl;

	FILE *file = fopen("res.txt", "w+");
	int *A = (int *)calloc(3, sizeof(int));
	int *m = (int *)calloc(6, sizeof(int));
	A[0] = 1;
	A[1] = 10;
	A[2] = 1000;
	m[0] = 1;
	m[1] = 1;
	m[2] = 2;
	m[3] = 2;
	m[4] = 2;
	m[5] = 1;

	for(int i = 0; i < 6; i++) {
		for(int j = 0; j < 3; j++) {
			fprintf(file, "%d ", i + 1);
			for(int k = 0; k < 4; k++) {
				if(E_n(i + 1, k + 1, A[j]) < 1e+10) {
					fprintf(file, "%lf ", E_n(i + 1, k + 1, A[j]));
					continue;
				}

				fprintf(file, "large ");
			}
			fprintf(file, "%d %d\n", m[i], A[j]);
		}
	}

	fclose(file);
	free(A);
	free(m);

	return 0;

}
