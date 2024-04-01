typedef double (* function_pointer)(double x);

double *Fourier_2d(function_pointer f, double p, int N);

double *error(double *Y, double *f, double p, int N);

double fourier(double *coeff, int num, double x);

double Lambda(int n, double p, int N);

double Scalar_(double *func, double h, int m, int n);

double *coeff_out(double *func, int values_number);
