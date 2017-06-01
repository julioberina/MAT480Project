#include <cmath>
#include <cstdio>
#define E 2.718281828
using namespace std;

// BEGIN FUNCTION 1 IMPLEMENTATIONS

double f1(double x, double y)
{
	return ((1.0 / 4.0)*pow(x, 4) - (1.0 / 3.0)*pow(x, 3) - 3 * pow(x, 2) + pow(y, 4) + 1.0);
}

// first component of gradient
double df1x(double x, double y)
{
	return (pow(x, 3) - pow(x, 2) - 6*x);
}

// second component of gradient
double df1y(double x, double y)
{
	return 4 * pow(y, 3);
}

// upper-left of Hessian
double df1xx(double x, double y)
{
	return (3 * pow(x, 2) - 2 * x - 6);
}

// upper-right of Hessian
double df1yx(double x, double y)
{
	return 0;
}

// lower-left of Hessian
double df1xy(double x, double y)
{
	return 0;
}

// lower-right of Hessian
double df1yy(double x, double y)
{
	return 12 * pow(y, 2);
}

// END FUNCTION 1 IMPLEMENTATIONS

// BEGIN FUNCTION 2 IMPLEMENTATIONS

double f2(double x, double y)
{
	return pow((-1)*E, ((-1)*pow(x, 2) - pow(y, 2)));
}

// first component of gradient
double df2x(double x, double y)
{
	return (2 * x*pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// second component of gradient
double df2y(double x, double y)
{
	return (2 * y*pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// upper-left of Hessian
double df2xx(double x, double y)
{
	return (2 * pow(E, ((-1)*pow(x, 2) - pow(y, 2))) - 4 * pow(x, 2)*pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// upper-right of Hessian
double df2yx(double x, double y)
{
	return (-4 * x * y * pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// lower-left of Hessian
double df2xy(double x, double y)
{
	return (-4 * x * y * pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// lower-right of Hessian
double df2yy(double x, double y)
{
	return (2 * pow(E, ((-1)*pow(x, 2) - pow(y, 2))) - 4 * pow(y, 2)*pow(E, ((-1)*pow(x, 2) - pow(y, 2))));
}

// END FUNCTION 2 IMPLEMENTATIONS



int main()
{

	return 0;
}