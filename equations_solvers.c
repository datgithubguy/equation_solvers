#include <math.h>
#include <stdio.h>

double		resolver_cardan(double p, double q);
double		solve_cardan(double a, double b, double c, double d);
double		solve_quartic(double a, double b, double c, double d, double e);

double		solve_cardan(double a, double b, double c, double d)
{
	double   p;
	double   q;
	double   ret;

	p = c / a - pow(b, 2) / (3.0 * pow(a, 2));
	q = (2.0 * pow(b, 3)) / (27.0 * pow(a, 3)) - (b * c) / (3.0 * pow(a, 2)) + d / a;
	ret = resolver_cardan(p, q);
	ret = ret -(b / (3.0 * a));
	return (ret);
}

double   resolver_cardan(double p, double q)
{
	double   root;
	double   delta;
	double   rho;
	double   theta;

	root = 0.0;
	delta = pow(q, 2) + (4.0 * pow(p, 3)) / 27.0;
	if (delta < 0)
	{
		rho = sqrt(pow(q, 2) - delta) / 2.0;
		theta = acos(-q / (2.0 * rho));
		root = 2.0 * cbrt(rho) * cos(theta / 3.0) - 2.0;
	}
	else if (delta == 0.0)
		root = cbrt(-4.0 * q);
	else
		root = cbrt((-q + sqrt(delta)) / 2.0) + cbrt((-q - sqrt(delta)) / 2.0);
	return (root);
}

double	solve_quartic(double a, double b, double c, double d, double e)
{
	double	y;
	double	p;
	double	q;
	double	r;
	double	del_o;
	double	del_b;
	double	x0;
	double	x1;
	double	x2;
	double	x3;

	p = -(3.0 * b * b) / 8.0 + c;
	q = (pow(b, 3) / 8.0) - (b * c / 2.0) + d;
	r = -3.0 * (pow(b / 4.0, 4)) + (c * pow(b, 2)) / 16.0 - (b * d) / 4.0 + e;
	y = solve_cardan(8.0, -4.0 * p, 8.0 * r, 4.0 * r * p - pow(q, 2));
	x0 = 0;
	x1 = 0;
	x2 = 0;
	x3 = 0;

	if ((-p + 2 * y) < 0)
	{
		return (0);
	}
	del_o = pow((-p + 2.0 * y), 1) - 4.0 * (y + (-q / (2.0 * sqrt(-p + 2.0 * y))));
	del_b = pow((-p + 2.0 * y), 1) - 4.0 * (y - (-q / (2.0 * sqrt(-p + 2.0 * y))));
	if (del_o < 0)
	{del_o = (double)del_o;
		// printf("no solution del_o\n");
	}
	else if (del_o)
	{
		x0 = (-(-p + 2.0 * y) + sqrt(del_o)) / (2.0 * a);
		x1 = (-(-p + 2.0 * y) - sqrt(del_o)) / (2.0 * a);
	}
	else
		x0 = -(-p + 2.0 * y);
	if (del_b < 0)
	{del_o = (double)del_o;
		// printf("no solution del_b\n");
	}
	else if (del_b)
	{
		x2 = (-(-p + 2.0 * y) + sqrt(del_b)) / (2.0 * a);
		x3 = (-(-p + 2.0 * y) - sqrt(del_b)) / (2.0 * a);
	}
	else
		x2 = -(-p + 2.0 * y);

	return (best_of_four(x0, x1, x2, x3));
}

int		best_of_four(int one, int two, int three, int four)
{
	int	best;

	best = 0;
	if (one > 0)
	{
		best = one;
	}
	if (two > best)
	{
		best = two;
	}
	if (three > best)
	{
		best = three;
	}
	if (four > best)
	{
		best = four;
	}
	return (best);
}
