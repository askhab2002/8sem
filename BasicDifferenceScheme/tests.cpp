#include <stdio.h>
#include <iostream>
#include "methods.hpp"

using namespace std;

int main(void) {

	double A = 1;
	int n = 1;
	int number_method = 1;

	cout << "   Введите значение параметра А" << endl;
	cin >> A;

	cout << "   Введите значения порядка дробления узлов n" << endl;
	cin >> n;

	cout << "   Введите номер разностной схемы" << endl;
	cin >> number_method;

	cout << E_n(number_method, n, A) << endl;

	return 0;

}
