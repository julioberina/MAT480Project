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

double* gradient1(double x, double y)
{
	double answer[] = { df1x(x, y), df1y(x, y) };
	return answer;
}

double** hessian1(double x, double y)
{
	double** answer = new double*[2];
	for (int i = 0; i < 2; ++i)
		answer[i] = new double[2];

	answer[0][0] = df1xx(x, y);
	answer[0][1] = df1yx(x, y);
	answer[1][0] = df1xy(x, y); 
	answer[1][1] = df1yy(x, y);

	return answer;
}

double** inv_hessian1(double x, double y)
{
	double** answer = new double*[2];
	for (int i = 0; i < 2; ++i)
		answer[i] = new double[2];

	double det = 1 / ((df1xx(x, y)*df1yy(x, y)) - (df1yx(x, y)*df1xy(x, y)));

	answer[0][0] = det * df1yy(x, y);
	answer[0][1] = (-1) * det * df1yx(x, y);
	answer[1][0] = (-1) * det * df1xy(x, y);
	answer[1][1] = det * df1xx(x, y);

	return answer;
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

double* gradient2(double x, double y)
{
	double answer[] = { df2x(x, y), df2y(x, y) };
	return answer;
}

double** hessian2(double x, double y)
{
	double** answer = new double*[2];
	for (int i = 0; i < 2; ++i)
		answer[i] = new double[2];

	answer[0][0] = df2xx(x, y);
	answer[0][1] = df2yx(x, y);
	answer[1][0] = df2xy(x, y);
	answer[1][1] = df2yy(x, y);

	return answer;
}

double** inv_hessian2(double x, double y)
{
	double** answer = new double*[2];
	for (int i = 0; i < 2; ++i)
		answer[i] = new double[2];

	double det = 1 / ((df2xx(x, y)*df2yy(x, y)) - (df2yx(x, y)*df2xy(x, y)));

	answer[0][0] = det * df2yy(x, y);
	answer[0][1] = (-1) * det * df2yx(x, y);
	answer[1][0] = (-1) * det * df2xy(x, y);
	answer[1][1] = det * df2xx(x, y);

	return answer;
}

// END FUNCTION 2 IMPLEMENTATIONS

void NewtonsMethod(int k, double g_x, double g_y, int fn, FILE* file)
{
	fprintf(file, "Newton's Method:\n\n");
	fprintf(file, "k\t\tx^(k)\n");

	// code

	fprintf(file, "\n");
}

void Project(int k, double guess_x, double guess_y, int function, int method)
{
	FILE *oFile;

	if (oFile = fopen("results.txt", "r"))
	{
		fclose(oFile);
		oFile = fopen("results.txt", "a");
	}
	else
	{
		fclose(oFile);
		oFile = fopen("results.txt", "w");
	}

	if (method == 1) // Newton's Method
		NewtonsMethod(k, guess_x, guess_y, function, oFile);

	// Close the file after operations are done
	fclose(oFile);
}

int main()
{

	return 0;
}