#ifndef SUPPORT_H
#define SUPPORT_H

#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include <iostream>

const double Pi = 3.14159265358979323846;

double sq(double x)
{
	return x*x;
}

struct Point
{
	double x, y;
};

std::tuple<double, double> linReg(const std::vector<Point>& points)
{
	double mean_x{ 0.0 };
	double mean_y{ 0.0 };

	for (const auto& p : points) {
		mean_x += p.x;
		mean_y += p.y;
	}
	mean_x /= points.size();
	mean_y /= points.size();

	double cov{ 0.0 }, var{ 0.0 };

	for (const auto& p : points) {
		cov += (p.x - mean_x) * (p.y - mean_y);
		var += (p.x - mean_x) * (p.x - mean_x);
	}
	double slope = cov / var;
	double intercept = mean_y - slope * mean_x;

	return std::make_tuple(slope, intercept);
}

std::tuple<double, double> linReg(const std::vector<double>& X, const std::vector<double>& Y)
{
	double mean_x{ 0.0 };
	double mean_y{ 0.0 };

	for (const auto& x : X) {
		mean_x += x;
	}
	mean_x /= X.size();

	for (const auto& y : Y) {
		mean_y += y;
	}
	mean_y /= Y.size();

	double cov{ 0.0 }, var{ 0.0 };

	for (int i = 0; i < X.size(); ++i)
	{
		cov += (X[i] - mean_x) * (Y[i] - mean_y);
		var += sq(X[i] - mean_x);
	}

	double slope = cov / var;
	double intercept = mean_y - slope * mean_x;

	return std::make_tuple(slope, intercept);
}

struct ParamExp
{
	double A, k;

	double operator()(double x) const
	{
		return A * std::exp(k * x);
	}
};

ParamExp fitExp(const std::vector<Point>& points)
{
	// ln of y's
	std::vector<Point> lnPoints = points;

	std::transform(points.begin(), points.end(), lnPoints.begin(), [](const Point& point) {
		return Point{ point.x, std::log(point.y) };
	});

	double slope, intercept;

	std::tie(slope, intercept) = linReg(lnPoints);
	return ParamExp{ std::exp(intercept), slope };
}

ParamExp fitExp(std::vector<Point>&& points)
{
	std::for_each(points.begin(), points.end(), [](Point& p) {
		p.y = std::log(p.y);
	});

	double slope, intercept;

	std::tie(slope, intercept) = linReg(points);
	return ParamExp{ std::exp(intercept), slope };
}

template <class Fun>
double integrate(Fun f, double start, double end, double coef = 0.75)
{
	if (end < start)
		throw std::runtime_error("integrate: end < start");
	if (end == start)
		return 0.0;
	// trapezoidal
	double I{ 0.0 };
	double stepLength{ (end - start)/10000.0 };
	stepLength *= coef;

	for (double x{ start }; x <= end; x += stepLength)
		I += (f(x) + f(x + stepLength)) * 0.5 * stepLength;

	return I;
}

double c10(double x, double a, double k_a)
{
	return (1 - x*((1 - k_a) / a));
}

double c20(double x, double a, double b, double k_a, double k_b)
{
	return k_a * std::exp(std::log(k_b / k_a) * (x - a) / (b - a));
}

double fusedFun(double x, double a, double b, double k_a, double k_b)
{
	if (x == 0.0)
		return 1.0;
	if (0.0 < x && x <= a)
		return c10(x, a, k_a);

	return c20(x, a, b, k_a, k_b);
}

using Matrix = std::vector < std::vector<double> >;
using Vector = std::vector < double >;

Matrix subMatrix(Matrix& m, int ir, int jr, int ic, int jc)
{
	Matrix subm(jr - ir + 1, Vector(jc - ic + 1));

	for (int r = 0; r < subm.size(); ++r)
	{
		for (int c = 0; c < subm[r].size(); ++c)
			subm[r][c] = m[r + ir][c + ic];
	}
	return subm;
}

std::ostream& operator<< (std::ostream& s, Matrix& m)
{
	for (int r = 0; r < m.size(); ++r)
	{
		for (int c = 0; c < m[r].size(); ++c)
			std::cout << /*std::setw(10) <<*/ m[r][c] << ' ';
		std::cout << '\n';
	}
	return s;
}


double a_func(double l, double delta, double a0,
	double l0, double mm, double K1c_max,
	double sigma, double alpha, double beta)
{
	double nom = sq(K1c_max) / (Pi * sq(sigma));
	nom *= (1 - delta);

	double denum = nom;
	nom -= l;
	denum -= l0;

	double frac = nom / denum;
	frac = 1 - std::pow(frac, beta);
	frac = std::pow(frac, 1.0 / alpha);
	frac = 1 + (mm - 1) * frac;
	return a0 * frac;
}

double c1(double x, double s1, double s2)
{
	return s1 * std::exp(s2 * x);
}

double der_c1(double x, double s1, double s2)
{
	return s1 * s2 * std::exp(s2 * x);
}

double der2_c1(double x, double s1, double s2)
{
	return s1 * sq(s2) * std::exp(s2 * x);
}

double F_torr(double x)
{
	return 1.0;
}



#endif // SUPPORT_H