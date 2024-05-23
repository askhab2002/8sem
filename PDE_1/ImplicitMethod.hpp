typedef double (* function_pointer_1)(double x);

typedef double (* function_pointer_2)(double x, double y);


void Sweep(double *y, double *a, double *b, double *c, double *f, int N);

double **SweepMethod(function_pointer_2 f, function_pointer_1 p, function_pointer_1 u_0, int T, int N);
