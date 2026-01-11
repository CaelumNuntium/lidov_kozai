double gamma_from_xy(double e, double i);
double beta(double e, double i, double g);
void ei2xy(double e, double i, double* x, double* y);
void xy2ei(double x, double y, double* e, double* i);
int lidov_kozai_stationary_point(int type, double nu2, double gamma, double sign_i, double* x0, double* y0, double* g0, double* x_dot, double* y_dot, double* g_dot, double* Omega_dot);
void lidov_kozai_eq(double t, double* x, double* params, double* res);
void three_body_eq(double t, double* x, double* params, double* res);
void lidov_kozai_analytical(int nt, double dt, int m, double nu2, double a0, double e0, double i0, double g0, double Omega0, double* a, double* x, double* y, double* g, double* Omega);
void lidov_kozai_numerical(int nt, double dt, int m, double Gm0, double Gm1, double Gm2, double a_, double e_, double i_, double g_, double Omega_, double p_long_, double M_0, double a0, double e0, double i0, double g0, double Omega0, double p_long0, double M0, double* a, double* e, double* incl, double* g, double* Omega, double* p_long, double* M);
