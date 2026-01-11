#include <stdlib.h>
#include <math.h>
#include "collo.h"
#include "legendre.h"
#include "orbits.h"
#include "lidov_kozai.h"

#define S 8

int sign(double);

int sign(double a)
{
    return a > 0.0 ? 1 : (a < 0.0 ? -1 : 0);
}

double gamma_from_xy(double x, double y)
{
    return (1.0 - x) * (1.0 - y);
}

double beta(double e, double i, double g)
{
    double s, sin_g;
    s = sin(i);
    sin_g = sin(g);
    return e * e * (2.0 - 5.0 * s * s * sin_g * sin_g);
}

void ei2xy(double e, double i, double* x, double* y)
{
    double s;
    *x = e * e;
    s = sin(i);
    *y = s * s;
}

void xy2ei(double x, double y, double* e, double* i)
{
    *e = sqrt(x);
    *i = asin(sqrt(y));
}

int lidov_kozai_stationary_point(int type, double nu2, double gamma, double sign_i, double* x0, double* y0, double* g0, double* x_dot, double* y_dot, double* g_dot, double* Omega_dot)
{
    switch(type)
    {
        case 1:
            *x0 = 0.0;
            *x_dot = 0.0;
            *y0 = 1.0 - gamma;
            *y_dot = 0.0;
            *g0 = 0.0;
            *g_dot = 0.0;
            *Omega_dot = -sign(sign_i) * nu2 * sqrt(gamma);
            break;
        case 2:
            *x0 = 1.0 - gamma;
            *x_dot = 0.0;
            *y0 = 0.0;
            *y_dot = 0.0;
            *Omega_dot = 0.0;
            *g0 = 0.0;
            *g_dot = nu2 * sqrt(gamma);
            break;
        case 3:
            *g0 = M_PI_2;
            *g_dot = 0.0;
            *x0 = 1.0 - sqrt(5.0 * gamma / 3.0);
            *x_dot = 0.0;
            *y0 = 1.0 - sqrt(0.6 * gamma);
            *y_dot = 0.0;
            *Omega_dot = -sign(sign_i) * nu2 * sqrt(0.6) * (1.0 + 4.0 * (*x0));
            break;
        default:
            return 1;
    }
    return 0;
}

void lidov_kozai_eq(double t, double* x, double* params, double* res) // x = [a, x, y, g, Omega]; params = [nu2, sign_i]
{
    double nu2, sign_i;
    nu2 = params[0];
    sign_i = params[1];
    res[0] = 0.0;
    res[1] = 10.0 * nu2 * x[1] * x[2] * sqrt(1.0 - x[1]) * sin(x[3]) * cos(x[3]);
    res[2] = -10.0 * nu2 * x[1] * x[2] * (1.0 - x[2]) / sqrt(1.0 - x[1]) * sin(x[3]) * cos(x[3]);
    res[3] = nu2 / sqrt(1.0 - x[1]) * (2.0 - 2.0 * x[1] + 5.0 * (x[1] - x[2]) * sin(x[3]) * sin(x[3]));
    res[4] = sign(sign_i) * nu2 / sqrt(1.0 - x[1]) * sqrt(1.0 - x[2]) * (1.0 - x[1] + 5.0 * x[1] * sin(x[3]) * sin(x[3]));
}

void three_body_eq(double t, double* x, double* params, double* res) // x = [x0, y0, z0, x1, y1, z1, x2, y2, z2, ...]; params = [Gm0, Gm1, Gm2]
{
    double r01_cube, r02_cube, r12_cube, r01_sq, r02_sq, r12_sq;
    res[0] = x[9];
    res[1] = x[10];
    res[2] = x[11];
    res[3] = x[12];
    res[4] = x[13];
    res[5] = x[14];
    res[6] = x[15];
    res[7] = x[16];
    res[8] = x[17];
    r01_sq = (x[3] - x[0]) * (x[3] - x[0]) + (x[4] - x[1]) * (x[4] - x[1]) + (x[5] - x[2]) * (x[5] - x[2]);
    r01_cube = r01_sq * sqrt(r01_sq);
    r02_sq = (x[6] - x[0]) * (x[6] - x[0]) + (x[7] - x[1]) * (x[7] - x[1]) + (x[8] - x[2]) * (x[8] - x[2]);
    r02_cube = r02_sq * sqrt(r02_sq);
    r12_sq = (x[6] - x[3]) * (x[6] - x[3]) + (x[7] - x[4]) * (x[7] - x[4]) + (x[8] - x[5]) * (x[8] - x[5]);
    r12_cube = r12_sq * sqrt(r12_sq);
    res[9] = params[1] / r01_cube * (x[3] - x[0]) + params[2] / r02_cube * (x[6] - x[0]);
    res[10] = params[1] / r01_cube * (x[4] - x[1]) + params[2] / r02_cube * (x[7] - x[1]);
    res[11] = params[1] / r01_cube * (x[5] - x[2]) + params[2] / r02_cube * (x[8] - x[2]);
    res[12] = params[0] / r01_cube * (x[0] - x[3]) + params[2] / r12_cube * (x[6] - x[3]);
    res[13] = params[0] / r01_cube * (x[1] - x[4]) + params[2] / r12_cube * (x[7] - x[4]);
    res[14] = params[0] / r01_cube * (x[2] - x[5]) + params[2] / r12_cube * (x[8] - x[5]);
    res[15] = params[0] / r02_cube * (x[0] - x[6]) + params[1] / r12_cube * (x[3] - x[6]);
    res[16] = params[0] / r02_cube * (x[1] - x[7]) + params[1] / r12_cube * (x[4] - x[7]);
    res[17] = params[0] / r02_cube * (x[2] - x[8]) + params[1] / r12_cube * (x[5] - x[8]);
}

void lidov_kozai_analytical(int nt, double dt, int m, double nu2, double a0, double e0, double i0, double g0, double Omega0, double* a, double* x, double* y, double* g, double* Omega)
{
    int i;
    double sign_i, x0, y0;
    double x_init[5], params[2];
    double* c;
    double* res;
    double* res_cur;
    ei2xy(e0, i0, &x0, &y0);
    sign_i = i0 > M_PI_2 ? -1.0 : 1.0;
    x_init[0] = a0;
    x_init[1] = x0;
    x_init[2] = y0;
    x_init[3] = g0;
    x_init[4] = Omega0;
    params[0] = nu2;
    params[1] = sign_i;
    c = (double*)malloc(S * sizeof(double));
    lobatto(S, c);
    res = (double*)malloc(5 * (nt + 1) * sizeof(double));
    collo(5, nt * m, m, dt / m, lidov_kozai_eq, x_init, params, S, c, res);
    for(i = 0; i <= nt; i++)
    {
        res_cur = res + i * 5;
        a[i] = res_cur[0];
        x[i] = res_cur[1];
        y[i] = res_cur[2];
        g[i] = res_cur[3];
        Omega[i] = res_cur[4];
    }
    free(res);
    free(c);
}

void lidov_kozai_numerical(int nt, double dt, int m, double Gm0, double Gm1, double Gm2, double a_, double e_, double i_, double g_, double Omega_, double p_long_, double M_0, double a0, double e0, double i0, double g0, double Omega0, double p_long0, double M0, double* a, double* e, double* incl, double* g, double* Omega, double* p_long, double* M)
{
    int i;
    double kappa_, kappa;
    double* c;
    double* res;
    double* res_cur;
    double x_init[18], params[3];
    vector3 r_init, v_init, rho0, u0, r0_init, v0_init, r1_init, v1_init, r0, v0, r1, v1, r, v;
    kappa_ = sqrt(Gm0 + Gm1 + Gm2);
    kappa = sqrt(Gm0 + Gm1);
    params[0] = Gm0;
    params[1] = Gm1;
    params[2] = Gm2;
    coords_from_orbital_elements(kappa_, a_, e_, i_, Omega_, g_, p_long_, M_0, &rho0, &u0);
    coords_from_orbital_elements(kappa, a0, e0, i0, Omega0, g0, p_long0, M0, &r_init, &v_init);
    r0_init = sum(rho0, cdot3(-Gm1 / (Gm0 + Gm1), r_init));
    r1_init = sum(rho0, cdot3(Gm0 / (Gm0 + Gm1), r_init));
    v0_init = sum(u0, cdot3(-Gm1 / (Gm0 + Gm1), v_init));
    v1_init = sum(u0, cdot3(Gm0 / (Gm0 + Gm1), v_init));
    vector_to_array(r0_init, x_init);
    vector_to_array(r1_init, x_init + 3);
    vector_to_array(new_vector(0.0, 0.0, 0.0), x_init + 6);
    vector_to_array(v0_init, x_init + 9);
    vector_to_array(v1_init, x_init + 12);
    vector_to_array(new_vector(0.0, 0.0, 0.0), x_init + 15);
    c = (double*)malloc(S * sizeof(double));
    lobatto(S, c);
    res = (double*)malloc(18 * (nt + 1) * sizeof(double));
    collo(18, nt * m, m, dt / m, three_body_eq, x_init, params, S, c, res);
    free(c);
    for(i = 0; i <= nt; i++)
    {
        res_cur = res + i * 18;
        r0 = from_array(res_cur);
        r1 = from_array(res_cur + 3);
        v0 = from_array(res_cur + 9);
        v1 = from_array(res_cur + 12);
        r = sum(r1, negative(r0));
        v = sum(v1, negative(v0));
        orbital_elements(kappa, r, v, a + i, e + i, incl + i, Omega + i, g + i, p_long + i, M + i);
    }
    free(res);
}
