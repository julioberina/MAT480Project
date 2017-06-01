#pragma warning (disable : 4996)
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

void NewtonsMethod(int k, double g_x, double g_y, int fn, FILE* file)
{
	fprintf(file, "Newton's Method (Function %d):\n\n", fn);
	fprintf(file, "k\t\tx^(i)\t\t\t\tf(x^(i))\n");

	double xi[] = { g_x, g_y };

	// allocate memory for inverse hessian and derivative of the function (gradient)
	double **ih = new double*[2];
	for (int i = 0; i < 2; ++i)
		ih[i] = new double[2];

	double *df = new double[2];
	double det = 0;

	// fill in the numbers

	if (fn == 1)
	{
		df[0] = df1x(xi[0], xi[1]);
		df[1] = df1y(xi[0], xi[1]);
		det = 1 / ((df1xx(xi[0], xi[1])*df1yy(xi[0], xi[1])) - (df1yx(xi[0], xi[1])*df1xy(xi[0], xi[1])));
		ih[0][0] = det * df1yy(xi[0], xi[1]);
		ih[0][1] = (-1) * det * df1yx(xi[0], xi[1]);
		ih[1][0] = (-1) * det * df1xy(xi[0], xi[1]);
		ih[1][1] = det * df1xx(xi[0], xi[1]);
	}
	else if (fn == 2)
	{
		df[0] = df2x(xi[0], xi[1]);
		df[1] = df2y(xi[0], xi[1]);
		det = 1 / ((df2xx(xi[0], xi[1])*df2yy(xi[0], xi[1])) - (df2yx(xi[0], xi[1])*df2xy(xi[0], xi[1])));
		ih[0][0] = det * df2yy(xi[0], xi[1]);
		ih[0][1] = (-1) * det * df2yx(xi[0], xi[1]);
		ih[1][0] = (-1) * det * df2xy(xi[0], xi[1]);
		ih[1][1] = det * df2xx(xi[0], xi[1]);
	}

	for (int i = 0; i < k; ++i)
	{
		fprintf(file, "%d\t\t", (i+1));
		fprintf(file, "(%lf, %lf)\t\t", xi[0], xi[1]);
		
		if (fn == 1)
			fprintf(file, "%lf\n", f1(xi[0], xi[1]));
		else if (fn == 2)
			fprintf(file, "%lf\n", f2(xi[0], xi[1]));

		double s[] = { (ih[0][0] * df[0] + ih[0][1] * df[1]), (ih[1][0] * df[0] + ih[1][1] * df[1]) };
		xi[0] = xi[0] - s[0];
		xi[1] = xi[1] - s[1];

		// recalculate gradient and hessian
		if (fn == 1)
		{
			df[0] = df1x(xi[0], xi[1]);
			df[1] = df1y(xi[0], xi[1]);
			det = 1 / ((df1xx(xi[0], xi[1])*df1yy(xi[0], xi[1])) - (df1yx(xi[0], xi[1])*df1xy(xi[0], xi[1])));
			ih[0][0] = det * df1yy(xi[0], xi[1]);
			ih[0][1] = (-1) * det * df1yx(xi[0], xi[1]);
			ih[1][0] = (-1) * det * df1xy(xi[0], xi[1]);
			ih[1][1] = det * df1xx(xi[0], xi[1]);
		}
		else if (fn == 2)
		{
			df[0] = df2x(xi[0], xi[1]);
			df[1] = df2y(xi[0], xi[1]);
			det = 1 / ((df2xx(xi[0], xi[1])*df2yy(xi[0], xi[1])) - (df2yx(xi[0], xi[1])*df2xy(xi[0], xi[1])));
			ih[0][0] = det * df2yy(xi[0], xi[1]);
			ih[0][1] = (-1) * det * df2yx(xi[0], xi[1]);
			ih[1][0] = (-1) * det * df2xy(xi[0], xi[1]);
			ih[1][1] = det * df2xx(xi[0], xi[1]);
		}
	}

	fprintf(file, "\n");

	delete[] df;
	for (int i = 0; i < 2; ++i)
		delete[] ih[i];
}

void SteepestDescent(int k, double g_x, double g_y, int fn, FILE* file)
{
	fprintf(file, "Steepest Descent (Function %d):\n\n", fn);
	fprintf(file, "k\t\tx^(i)\t\t\t\tf(x^(i))\n");

	double* df = new double[2];
	double xi[] = { g_x, g_y };
	double magnitude = 0.0;

	// calculate gradient based on function
	if (fn == 1)
	{
		df[0] = df1x(xi[0], xi[1]);
		df[1] = df1y(xi[0], xi[1]);
	}
	else if (fn == 2)
	{
		df[0] = df2x(xi[0], xi[1]);
		df[1] = df2y(xi[0], xi[1]);
	}

	for (int i = 0; i < k; ++i)
	{
		fprintf(file, "%d\t\t", (i + 1));
		fprintf(file, "(%lf, %lf)\t\t", xi[0], xi[1]);

		if (fn == 1)
			fprintf(file, "%lf\n", f1(xi[0], xi[1]));
		else if (fn == 2)
			fprintf(file, "%lf\n", f2(xi[0], xi[1]));

		xi[0] = xi[0] - 0.1*df[0];
		xi[1] = xi[1] - 0.1*df[1];

		magnitude = sqrt(pow(df[0], 2) + pow(df[1], 2));
		if (magnitude <= 0.1)
			break;
		else
		{
			// recalculate gradient based on function
			if (fn == 1)
			{
				df[0] = df1x(xi[0], xi[1]);
				df[1] = df1y(xi[0], xi[1]);
			}
			else if (fn == 2)
			{
				df[0] = df2x(xi[0], xi[1]);
				df[1] = df2y(xi[0], xi[1]);
			}
		}
	}

	fprintf(file, "\n");

	delete[] df;
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
		oFile = fopen("results.txt", "w");

	switch (method)
	{
	case 1: // Newton's Method
		NewtonsMethod(k, guess_x, guess_y, function, oFile);
		break;
	case 2: // Steepest Descent
		SteepestDescent(k, guess_x, guess_y, function, oFile);
		break;
	}

	// Close the file after operations are done
	fclose(oFile);
}

int main()
{
	Project(20, 2.0, 1.0, 1, 1);
	Project(20, 2.0, 1.0, 1, 2);

	Project(20, 0.2, 0.4, 2, 1);
	Project(20, 0.2, 0.4, 2, 2);

	// Own guesses
	// Project(20, 1.0, 3.0, 1, 1); // [1, 3]
	// Project(20, 0.1, 0.1, 2, 1); // [0.8, 0.9]

	return 0;
}