typedef double (* function_pointer_2d)(double x, double y);

double **coeff_out_2d(int node_number, function_pointer_2d func);

double fourier(double **coeff, int num, double x, double y);
