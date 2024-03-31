typedef void  (* function_pointer)(double x, double *y, int n, double *temporary);

void K(function_pointer f, double x, double *y, double h, int n, double *k);

void Yn(function_pointer f, double x, double *y, double h, int n, double *y_n, double *k1, double *k2, double *k3, double *k4);

void En(function_pointer f, double x, double *y, double h, int n, double *E_n, double *k1, double *k2, double *k3, double *k4);

double L2_h(double *x, double *y, double h, int n);

double L0(double *x, int n);

double L0(double *x, double *y, int n);

double L2_h_diff(double *x, double *y, double h, int n);
