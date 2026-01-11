#include <math.h>
#include "vector.h"
#include "orbits.h"

#define MIN_ECCENTRICITY 1e-10
#define MIN_INCLINATION 1e-10

void coords_from_orbital_elements(double kappa, double a, double e, double i, double Omega, double omega, double p_long, double M, vector3* r, vector3* v)
{
	const double eps = 2.5e-16;
	const int NMAX = 1000;
	int k;
	double E, E_last, n;
	double KAPPASQ = kappa * kappa;
	vector3 r0, v0;
	matrix33 S;
	E_last = M;
	E = M + e * sin(M);
	k = 0;
	while (fabs(E - E_last) > eps && k < NMAX)
	{
		E_last = E;
		E = M + e * sin(E_last);
		k++;
		if(k > NMAX)
		{
			break;
		}
	}
	n = sqrt(KAPPASQ) * pow(a, -1.5);
	r0.x = a * (cos(E) - e);
	r0.y = a * sqrt(1 - e * e) * sin(E);
	r0.z = 0.0;
	v0.x = -n * a * a / norm3(r0) * sin(E);
	v0.y = n * a * a / norm3(r0) * sqrt(1.0 - e * e) * cos(E);
	v0.z = 0.0;
	if(fabs(i) < MIN_INCLINATION && e < MIN_ECCENTRICITY)
	{
		*r = r0;
		*v = v0;
		return;
	}
	else if(fabs(i) < MIN_INCLINATION)
	{
		S = new_matrix(cos(p_long), -sin(p_long), 0.0, \
			sin(p_long), cos(p_long), 0.0, \
			0.0, 0.0, 1.0);
	}
	else if(e < MIN_ECCENTRICITY)
	{
		S = new_matrix(cos(Omega), -sin(Omega) * cos(i), sin(Omega) * sin(i), \
			sin(Omega), cos(Omega) * cos(i), -cos(Omega) * sin(i), \
			0.0, sin(i), cos(i));
	}
	else
	{
		S = new_matrix(\
			cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i), \
			-cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i), \
			sin(Omega) * sin(i), \
			sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i), \
			-sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i), \
			-cos(Omega) * sin(i), \
			sin(omega) * sin(i), \
			cos(omega) * sin(i), \
			cos(i));
	}
	*r = linmap(S, r0);
	*v = linmap(S, v0);
}

void orbital_elements(double kappa, vector3 r, vector3 v, double* a, double* e, double* i, double* Omega, double* omega, double* p_long, double* M)
{
	double KAPPASQ = kappa * kappa;
	vector3 n;
	double c, p, rr, u, theta, dr, E, beta;
	rr = norm3(r);
	*a = 1.0 / (2.0 / rr - dot(v, v) / KAPPASQ);
	n = cross(r, v);
	c = norm3(n);
	n = cdot3(1.0 / c, n);
	*i = atan2(sqrt(n.x * n.x + n.y * n.y), n.z);
	p = c * c / KAPPASQ;
	if(p / *a < 1.0)
	{
		*e = sqrt(1.0 - p / *a);
	}
	else
	{
		*e = 0.0;
	}
	dr = dot(r, v) / rr;
	if(fabs(*i) < MIN_INCLINATION && *e < MIN_ECCENTRICITY)
	{
		*Omega = 0.0;
		*omega = 0.0;
		*M = to_interval(atan2(r.y, r.x), 0.0, 2.0 * M_PI);
		*p_long = 0.0;
	}
	else if(fabs(*i) < MIN_INCLINATION)
	{
		*Omega = 0.0;
		theta = atan2(dr * sqrt(p) / kappa, p / rr - 1.0);
		E = to_interval(2.0 * atan(sqrt((1.0 - *e) / (1.0 + *e)) * tan(0.5 * theta)), 0.0, 2.0 * M_PI);
		*M = to_interval(E - *e * sin(E), 0.0, 2.0 * M_PI);
		*p_long = to_interval(atan2(r.y, r.x) - theta, 0.0, 2.0 * M_PI);
		*omega = *p_long;
	}
	else if(*e < MIN_ECCENTRICITY)
	{
		*Omega = to_interval(atan2(n.x, -n.y), 0.0, 2.0 * M_PI);
		*omega = 0.0;
		u = atan2(r.z / sin(*i), r.x * cos(*Omega) + r.y * sin(*Omega));
		*M = to_interval(u, 0.0, 2.0 * M_PI);
		*p_long = *Omega;
	}
	else
	{
		*Omega = to_interval(atan2(n.x, -n.y), 0.0, 2.0 * M_PI);
		u = atan2(r.z / sin(*i), r.x * cos(*Omega) + r.y * sin(*Omega));
		theta = atan2(dr * sqrt(p) / kappa, p / rr - 1.0);
		*omega = to_interval(u - theta, 0.0, 2.0 * M_PI);
		E = to_interval(2.0 * atan(sqrt((1.0 - *e) / (1.0 + *e)) * tan(0.5 * theta)), 0.0, 2.0 * M_PI);
		*M = to_interval(E - *e * sin(E), 0.0, 2.0 * M_PI);
		*p_long = to_interval(*Omega + *omega, 0.0, 2.0 * M_PI);
	}
	//beta = *e / (1.0 - (*e) * (*e));
	//E = theta - 2.0 * atan(beta * sin(theta) / (1.0 + beta * cos(theta)));
}

double degrees(double radians)
{
    return radians * 180.0 / M_PI;
}

double radians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double to_interval(double value, double left_bound, double right_bound)
{
	double len, periods, res;
    len = right_bound - left_bound;
    periods = (value - left_bound) / len;
	res = value - floor(periods) * len;
	if(res < left_bound)
	{
		res += len;
	}
	if(res >= right_bound)
	{
		res -= len;
	}
    return res;
}
