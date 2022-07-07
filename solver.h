void gen_Z_drain(int nx, double * xs, double * Z, double height, double centre);

void apply_topography(double * Z_arr, double * xs, double k, double gamma_of_k, double theta, int nx);

void WB_solver(double ** W, double ** W_forward, double ** W_L_stars, double ** W_R_stars, double * lmbda_neg, double * lmbda_pos, double * Z, int nx, double dx, double T_stop, double * xs, double k, double gamma_of_k, double sig_theta, FILE * CURVE_DATA, FILE * TOP_DATA, FILE * TIMES);