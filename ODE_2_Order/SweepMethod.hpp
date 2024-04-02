typedef double (* function_pointer)(double x);

double *SweepMethod(function_pointer f, function_pointer p, int N);

double *Sweep(double *a, double *b, double *c, double *f, int N);
