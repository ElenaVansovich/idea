#include <iostream>
#include "math.h"

using namespace std;

double f(double x)
{
	return x*x*x - 2048.125*x*x + 1048832 * x - 131072;
}

double simpleIter(double x, double e)
{
	double result = 0;
	while (true)
	{
		result = x;
		x = (-1)*((x*x*x - 2048.125*x*x - 131072) / 1048832);
		cout << "е " << abs(x - result) << endl;
		if (abs(x - result) <= e)
		{
			break;
		}

	}
	return x;
}

double secant(double x, double e, double q)
{
	double f1;
	double result = 0;

	while (true)
	{
		result = x;
		f1 = (f(x) - f(x - q)) / q;
		x = x - f(x) / (f1);
		cout << "е " << abs(x - result) << endl;
		if (abs(x - result) <= e)
		{
			break;
		}

	}
	return x;
}

void matrJacobi(double* x, double** result)
{
	result[0][0] = 1;
	result[0][1] = 10 * x[1] - 3 * x[1] * x[1] - 2;
	result[1][0] = 1;
	result[1][1] = 2 * x[1] + 3 * x[1] * x[1] - 14;
}

void vectorFunc(double* x, double* result)
{
	result[0] = (-1)*(x[0] + x[1] * (x[1] * (5 - x[1]) - 2) - 13);
	result[1] = (-1)*(x[0] + x[1] * (x[1] * (1 + x[1]) - 14) - 29);
}

void inversMatr(double** x, double** result)
{
	double tmpDet = (1 / abs(x[0][0] * x[1][1] - x[0][1] * x[1][0]));
	
	result[0][0] = tmpDet * x[1][1];
	result[0][1] = tmpDet * x[0][1] * (-1);
	result[1][0] = tmpDet * x[1][0] * (-1);
	result[1][1] = tmpDet * x[0][0];
}

void funcDX(double** wx, double* fx, double* result)//поправка
{
	double* x = new double[2];

	double** mass = new double*[2];
	for (int i = 0; i < 2; i++)
	{
		mass[i] = new double[2];
	}
	inversMatr(wx, mass);

	result[0] = mass[0][0] * fx[0] + mass[0][1] * fx[1];
	result[1] = mass[1][0] * fx[0] + mass[1][1] * fx[1];

	delete[] x;
	for (int i = 0; i < 2; i++){
		delete[] mass[i];
	}
	delete[] mass;
}

void vectorSum(double* A, double* B)
{
	A[0] += B[0];
	A[1] += B[1];
}

double maxNorm(double* x)
{
	if (abs(x[0]) > abs(x[1]))
		return abs(x[0]);
	else
		return abs(x[1]);
}

void newton(double* x, double e)
{
	double** wх = new double*[2]{};
	for (int i = 0; i < 2; i++)
	{
		wх[i] = new double[2]{};
	}

	double* fx = new double[2]{};
	double* dx = new double[2]{};
	while (true)
	{
		matrJacobi(x, wх);

		vectorFunc(x, fx);

		funcDX(wх, fx, dx);

		vectorSum(x, dx);

		cout << "е " << maxNorm(dx) << endl;
		if (maxNorm(dx) <= e)
		{
			break;
		}
	}

	for (int i = 0; i < 2; i++){
		delete[] wх[i];
	}
	delete[] wх;
	delete[] fx;
	delete[] dx;
}

int main()
{
	setlocale(LC_ALL, "rus");
	double x = 0;
	double result = 0;
	double e = pow(10, -10);

	cout << "Метод простой итерации:" << endl;
	cout << "Первый корень: " << simpleIter(x, e) << endl;
	x = 1024;
	cout << "Второй корень: " << simpleIter(x, e) << endl;

	double q = 0.1;
	x = 0;
	cout << "Метод секущих:" << endl;
	cout << "Первый корень: "  << secant(x, e, q) << endl;
	x = 1024;
	cout << "Второй корень: " << secant(x, e, q) << endl;
	
	cout << "Метод Ньютона для системы:" << endl;
	double* mass = new double[2]{};
	mass[0] = 15;
	mass[1] = -2;
	newton(mass, pow(10, -10));
	cout << mass[0] << " " << mass[1] << endl;
	delete[] mass;
	return 0;
}