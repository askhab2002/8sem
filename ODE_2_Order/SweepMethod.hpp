typedef double (* function_pointer)(double x);

double *SweepMethod(function_pointer p, function_pointer f, int N);

double *Sweep(double *a, double *b, double *c, double *f, int N);
