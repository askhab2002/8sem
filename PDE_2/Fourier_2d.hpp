typedef double (* function_pointer_2d)(double x, double y);

void coeff_out_2d(int node_number, double **func, double **coeff_new);

double fourier(double **coeff, int num, double x, double y);

double Lambda_out(int n, double p, int N);

double *generate_equel(double left, int num);

double **do_matrix(function_pointer_2d func, double *nodes, int node_number);
