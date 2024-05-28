
typedef double (* function_pointer)(double x);

double **DoMatrix(int N, double p);
double *DoB(function_pointer f, int N);
double *DoP(int N);
double **DoMatrixP(int N, double *p);
double *DoF(function_pointer f, int N);


double Tau(double **A, int N, double *q);
double *StepB(double **A, double **inv_A, double *b, double *x_k, int N, double Tau);
double BSolver(double *x, double **A, double *b, double Tau, int N, double eps, int mIter);


double **MatrixMult(double **A, double **B, int N);
double *MatrixVector(double **A, double *x, int N);
double **Inverse(int N);


double fourier(double *coeff, int num, double x);
double Lambda(int n, double p, int N);
double Lambda_k(int n, int k, int N, double h);
double cosinus(double x);
double cosi(int m, int k, int N);
double Scalar_(double *func, function_pointer cosin, double h, int m, int n);
double *coeff_out(double *func, function_pointer cosin, int values_number);


void test_3(int N, function_pointer f, int mIter);
double error_B_(double *x, double *b, double **A, int N);
double *TriangularSystemP(int N, double *b, function_pointer cosin);
double *error_B(double *Y, double *f, double h, int N);

