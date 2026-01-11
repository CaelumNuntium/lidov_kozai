#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "orbits.h"
#include "lidov_kozai.h"

#define N_PARAMS 16

int write_result(const char*, int, double, int, double*, double*, double*, double*, double*, double*);

int main()
{
    const char* params[N_PARAMS] = { "stat_point_type", "N_steps", "m_steps", "delta_t", "Gm0", "Gm1", "Gm2", "a'", "e'", "i'", "g'", "Omega'", "M'(0)", "a(0)", "M(0)", "gamma" };
    char* values[N_PARAMS];
    int i, stat_point_type, n_steps, m;
    double nu, nu1, nu2, gamma, sign_i, x0, y0, g0, Omega0, x_dot, y_dot, g_dot, Omega_dot, dt, tt;
    double Gm0, Gm1, Gm2, a_, e_, i_, g_, Omega_, M0_, a0, e0, i0, M0;
    double* x_res_a;
    double* y_res_a;
    double* g_res_a;
    double* Omega_res_a;
    double* a_res_a;
    double* p_long_res_a;
    double* x_res_n;
    double* y_res_n;
    double* g_res_n;
    double* Omega_res_n;
    double* a_res_n;
    double* e_res_n;
    double* i_res_n;
    double* M_res_n;
    double* p_long_res_n;
    for(i = 0; i < N_PARAMS; i++)
    {
        values[i] = (char*)malloc(100 * sizeof(char));
    }
    read_config("config.ini", N_PARAMS, params, values, NULL);
    sscanf(values[0], "%d", &stat_point_type);
    sscanf(values[1], "%d", &n_steps);
    sscanf(values[2], "%d", &m);
    sscanf(values[3], "%lf", &dt);
    sscanf(values[4], "%lf", &Gm0);
    sscanf(values[5], "%lf", &Gm1);
    sscanf(values[6], "%lf", &Gm2);
    sscanf(values[7], "%lf", &a_);
    sscanf(values[8], "%lf", &e_);
    sscanf(values[9], "%lf", &i_);
    sscanf(values[10], "%lf", &g_);
    sscanf(values[11], "%lf", &Omega_);
    sscanf(values[12], "%lf", &M0_);
    sscanf(values[13], "%lf", &a0);
    sscanf(values[14], "%lf", &M0);
    sscanf(values[15], "%lf", &gamma);
    sign_i = 1;
    printf("N_steps = %d\nm_steps = %d\ndelta_t = %lf\n\n", n_steps, m, dt);
    printf("stat_point_type = %d\ngamma = %lf\n", stat_point_type, gamma);
    printf("Gm0 = %lf\nGm1 = %lf\nGm2 = %lf\n\n", Gm0, Gm1, Gm2);
    printf("a' = %lf\ne' = %lf\ni' = %lf\ng' = %lf\nOmega' = %lf\nM'(0) = %lf\n\n", a_, e_, i_, g_, Omega_, M0_);
    printf("a(0) = %lf\nM(0) = %lf\n\n", a0, M0);
    nu = Gm2 / (Gm0 + Gm1);
    nu1 = nu * a0 * a0 * a0 / (a_ * a_ * a_) * 0.375 / ((1.0 - e_ * e_) * sqrt(1.0 - e_ * e_));
    nu2 = 2.0 * sqrt(Gm0 + Gm1) * nu1 / (a0 * sqrt(a0));
    printf("nu2 = %lf\n", nu2); 
    for(i = 0; i < N_PARAMS; i++)
    {
        free(values[i]);
    }
    if(stat_point_type < 1 || stat_point_type > 3)
    {
        fprintf(stderr, "Error: Invalid value: stat_point_type = %d\n", stat_point_type);
        return 1;
    }
    lidov_kozai_stationary_point(stat_point_type, nu2, gamma, sign_i, &x0, &y0, &g0, &x_dot, &y_dot, &g_dot, &Omega_dot);
    printf("x(0) = %lf\ny(0) = %lf\ng(0) = %lf\n", x0, y0, g0);
    x_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    y_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    g_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    Omega_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    a_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    p_long_res_a = (double*)malloc((n_steps + 1) * sizeof(double));
    for(i = 0; i <= n_steps; i++)
    {
        tt = dt * i;
        a_res_a[i] = a0;
        x_res_a[i] = x0 + x_dot * tt;
        y_res_a[i] = y0 + y_dot * tt;
        g_res_a[i] = g0 + g_dot * tt;
        Omega_res_a[i] = Omega0 + Omega_dot * tt;
        p_long_res_a[i] = g_res_a[i] + Omega_res_a[i];
    }
    write_result("result_analytical.dat", n_steps, dt, stat_point_type, a_res_a, x_res_a, y_res_a, g_res_a, Omega_res_a, p_long_res_a);
    free(a_res_a);
    free(x_res_a);
    free(y_res_a);
    free(g_res_a);
    free(Omega_res_a);
    free(p_long_res_a);
    xy2ei(x0, y0, &e0, &i0);
    printf("e(0) = %lf\ni(0) = %lf\n", e0, i0);
    a_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    e_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    i_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    g_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    Omega_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    M_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    p_long_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    lidov_kozai_numerical(n_steps, dt, m, Gm0, Gm1, Gm2, a_, e_, i_, g_, Omega_, g_ + Omega_, M0_, a0, e0, i0, g0, Omega0, g0 + Omega0, M0, a_res_n, e_res_n, i_res_n, g_res_n, Omega_res_n, p_long_res_n, M_res_n);
    free(M_res_n);
    x_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    y_res_n = (double*)malloc((n_steps + 1) * sizeof(double));
    for(i = 0; i <= n_steps; i++)
    {
        ei2xy(e_res_n[i], i_res_n[i], x_res_n + i, y_res_n + i);
    }
    free(e_res_n);
    free(i_res_n);
    write_result("result_numerical.dat", n_steps, dt, stat_point_type, a_res_n, x_res_n, y_res_n, g_res_n, Omega_res_n, p_long_res_n); 
    free(a_res_n);
    free(x_res_n);
    free(y_res_n);
    free(g_res_n);
    free(Omega_res_n);
    free(p_long_res_n);
    return 0;
}

int write_result(const char* filename, int n, double dt, int stat_point_type, double* a, double* x, double* y, double* g, double* Omega, double* p_long)
{
    int i;
    FILE* out;
    if(!(out = fopen(filename, "w")))
    {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return 1;
    }
    fprintf(out, "#     t            a         x=e^2     y=(sin(i))^2     g,deg     Omega,deg       pericenter_longitude\n");
    for(i = 0; i <= n; i++)
    {
        fprintf(out, "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", i * dt, a[i], x[i], y[i], to_interval(degrees(g[i]), 0.0, 360.0), to_interval(degrees(Omega[i]), 0.0, 360.0), to_interval(degrees(p_long[i]), 0.0, 360.0));
    }
    fclose(out);
    return 0;
}
