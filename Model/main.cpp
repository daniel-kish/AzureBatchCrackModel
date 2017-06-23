#include "Support.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <chrono>
#include <string>
#include <chrono>
#include <sstream>

using namespace std;
using namespace chrono;

// global parameters
double D{ 10E-11 };
double V_H{ 1.96E-6 };
double R{ 8.314 };
double T{ 293.0 };
double C{ 1.095E-12 };
double n{ 3.24 };
double w{ 50E-3 };


double Speed_fatigue(double K_1, double ff)
{
	return C * ff * std::pow(K_1, n);
}

double K1_(double sigma, double l)
{
	return sigma * std::sqrt(Pi * l) /** F_torr(l)*/;
}


//
//double c10(double x, double a, double k_a)
//{
//	return (1 - x*((1 - k_a) / a));
//}
//
//double c20(double x, double a, double b, double k_a, double k_b)
//{
//	return k_a * std::exp(std::log(k_b / k_a) * (x - a) / (b - a));
//}

double f(double x, double a, double b, double k_a, double k_b)
{
	if (x == 0.0)
		return 1.0;
	if (0.0 < x && x <= a)
		return c10(x, a, k_a);

	return c20(x, a, b, k_a, k_b);
}




Matrix Solve_env(double l0, double delta, double a0, double sigma, double K_1cmin,
	double K1c_max, double alpha, double beta, double mm, double alpha_2, double beta_2, double omega)
{
	int sz = 5000;
	Vector l(sz); l[0] = l0;
	Vector k_a(sz); k_a[0] = 0.0025;
	Vector k_b(sz); k_b[0] = 0.000001;
	Matrix aaa(sz, Vector(20)); // Matrix sz x sz
	int i = 0;
	int j = 0;
	Vector a(sz);
	Vector b(sz);
	Vector SS_1(sz), SS_2(sz);
	double m{ 0.0 };
	double B{ 0.0 };
	double deltaT{ 0.0 };

	auto cond = [&](double l_i){
		return l_i < ((1 - delta) * (sq(K1c_max)) / (Pi*sq(sigma)));
	};
	auto expParams = [&](int i){
		int nr_points = int(5.0 * b[i] / a[i]);
		double step = a[i] * 0.2;

		Vector x(nr_points);

		for (int i_i = 1; i_i < nr_points; ++i_i)
			x[i_i] = x[i_i - 1] + step;

		Vector F_vec(x.size());

		std::transform(std::begin(x), std::end(x), std::begin(F_vec), [&](double& elem) {
			return std::log(f(elem, a[i], b[i], k_a[i], k_b[i]));
		});

		double slope, intercept;
		std::tie(slope, intercept) = linReg(x, F_vec);
		return std::make_tuple(std::exp(intercept), slope);
	};
	int x;
	auto m_formula = [&] {
		double num{ 0.0 };
		num = (K1_(sigma, l[i]) - K_1cmin) / (K1c_max - K_1cmin);
		num = std::pow(num, alpha);
		num = 1.0 - num;
		num = std::pow(num, 1.0 / beta);
		num /= omega;

		double denom = integrate(
			[&](double x) { return c1(x, SS_1[i], SS_2[i]); }, 0.0, a[i]
			);
		denom /= a[i];
		return num / denom;
	};
	auto B_formula = [&] {
		auto integrand = [&](double x){
			double tmp{ 0.0 };
			tmp = der_c1(x, SS_1[i], SS_2[i]);
			tmp *= c1(x, SS_1[i], SS_2[i]);
			tmp /= std::sqrt(Pi) * 2.0 * x * std::sqrt(x);
			return tmp;
		};
		return integrate(integrand, a[i], b[i]);
	};
	auto deltaT_formula = [&]{
		auto num = [&]{
			return std::log(m) * integrate([&](double x){ return sq(c1(x, SS_1[i], SS_2[i])); },
				a[i], b[i], 0.4);
		};
		double Num = num();
		//double tmp1 = num();
		auto denom = [&]{
			double tmp = K1_(sigma, l[i]);
			double first = sqrt(2.0 / Pi)*(-D * V_H) / (3.0 * R * T) * tmp * B;
			double second = D *
				integrate([&](double x) {
				return c1(x, SS_1[i], SS_2[i]) * der2_c1(x, SS_1[i], SS_2[i]);
			}, a[i], b[i], 0.4);
			return first + second;
		};
		double Denom = denom();
		//double tmp2 = denom();
		return num() / denom();
	};

	while (cond(l[i]))
	{
		a[i] = a_func(l[i], delta, a0, l0, mm, K1c_max, sigma, alpha_2, beta_2);
		b[i] = 500 * a[i];
		double expIntercept, slope;
		std::tie(expIntercept, slope) = expParams(i);
		SS_1[i] = expIntercept;
		SS_2[i] = slope;
		m = m_formula();
		B = B_formula();
		deltaT = deltaT_formula();

		aaa[i][6] = deltaT;
		aaa[i][11] = a[i] / aaa[i][6];

		aaa[i][1] = l[i];
		aaa[i][2] = K1_(sigma, l[i]);
		aaa[i][3] = m;
		aaa[i][4] = k_a[i];
		aaa[i][5] = k_b[i];

		if (i == 0)
			aaa[i + 1][7] = aaa[i][6];
		else
			aaa[i + 1][7] = aaa[i][7] + deltaT;

		aaa[i][8] = integrate([&](double x){ return c1(x, SS_1[i], SS_2[i]); }, 0.0, a[i]) / a[i];
		l[i + 1] = l[i] + a[i];
		a[i + 1] = a_func(l[i + 1], delta, a0, l0, mm, K1c_max, sigma, alpha_2, beta_2);
		b[i + 1] = 500.0 * a[i + 1];
		k_b[i + 1] = m * c1(a[i] + b[i + 1], SS_1[i], SS_2[i]);
		k_a[i + 1] = m * c1(a[i + 1] + a[i], SS_1[i], SS_2[i]);

		double expIntercept2, slope2;
		std::tie(expIntercept2, slope2) = expParams(i + 1);
		aaa[i][9] = a[i];
		aaa[i][10] = 0.0;
		aaa[i][12] = j;
		aaa[i][13] = SS_1[i];
		aaa[i][14] = SS_2[i];
		auto cond1 = [&](){
			double I = integrate([&](double x){ return c1(x, expIntercept2, expIntercept2);	}, 0.0, a[i + 1]) / a[i + 1];
			auto smthElse = [&](){
				double num{ 0.0 };
				num = (K1_(sigma, l[i + 1]) - K_1cmin) / (K1c_max - K_1cmin);
				num = std::pow(num, alpha);
				num = 1.0 - num;
				num = std::pow(num, 1.0 / beta);
				num /= omega;
				return num;
			};
			auto SmthElse = smthElse();
			return (deltaT <= 0.0 || I >= SmthElse);
		};
		if (cond1())
			break;
		i++;
	}
	if (i == 0) {
		aaa[0][14] = 0.0;
		aaa[0][1] = l0;
		return aaa;
	}
	else {
		aaa.resize(i);
		return aaa;
	}
}


Matrix Solve_cycle(double l0, double delta, double a0, double sigma, double Kmin, double Kmax,
	double alpha, double beta, double mm, double alpha2, double beta2, double fff, double omega)
{
	int sz{ 1400 };
	Vector l(sz);
	Vector ka(sz);
	Vector kb(sz);
	Vector a(sz);
	Vector b(sz);
	Vector expInter(sz);
	Vector slope(sz);
	double m{ 0.0 };
	double B{ 0.0 };
	double deltaT{ 0.0 };
	Matrix data(sz, Vector(20));

	int i = 1;
	l[i] = l0;
	ka[i] = 0.0025;
	kb[i] = 0.000001;
	int j = 0;
	int jjjprerup = 1;
	int kk = 0;
	int jjj{ 0 };

	auto Lmax = [&] {
		return (1.0 - delta) * sq(Kmax / sigma) / Pi;
	};

	auto expParams = [&](int i){
		int nr_points = int(5.0 * b[i] / a[i]);
		double step = a[i] * 0.2;

		Vector x(nr_points);

		for (int i_i = 1; i_i < nr_points; ++i_i)
			x[i_i] = x[i_i - 1] + step;

		Vector F_vec(x.size());

		std::transform(std::begin(x), std::end(x), std::begin(F_vec), [&](double& elem) {
			return std::log(f(elem, a[i], b[i], ka[i], kb[i]));
		});

		double slope, intercept;
		std::tie(slope, intercept) = linReg(x, F_vec);
		return std::make_tuple(std::exp(intercept), slope);
	};
	auto m_formula = [&] {
		double num{ 0.0 };
		num = (K1_(sigma, l[i]) - Kmin) / (Kmax - Kmin);
		num = std::pow(num, alpha);
		num = 1.0 - num;
		num = std::pow(num, 1.0 / beta);
		num /= omega;

		double denom = integrate(
			[&](double x) { return c1(x, expInter[i], slope[i]); }, 0.0, a[i]
			);
		denom /= a[i];
		return num / denom;
	};
	auto B_formula = [&] {
		auto integrand = [&](double x){
			double tmp{ 0.0 };
			tmp = der_c1(x, expInter[i], slope[i]);
			tmp *= c1(x, expInter[i], slope[i]);
			tmp /= std::sqrt(Pi) * 2.0 * x * std::sqrt(x);
			return tmp;
		};
		return integrate(integrand, a[i], b[i]);
	};
	auto deltaT_formula = [&]{
		auto num = [&]{
			return std::log(m) * integrate([&](double x){ return sq(c1(x, expInter[i], slope[i])); },
				a[i], b[i], 0.4);
		};
		double Num = num();
		//double tmp1 = num();
		auto denom = [&]{
			double tmp = K1_(sigma, l[i]);
			double first = 1.0E6 * (-D * V_H) / (3.0 * R * T) * tmp * B;
			double second = D *
				integrate([&](double x) {
				return c1(x, expInter[i], slope[i]) * der2_c1(x, expInter[i], slope[i]);
			}, a[i], b[i], 0.4);
			return first + second;
		};
		double Denom = denom();
		//double tmp2 = denom();
		return num() / denom();
	};
	auto m_formula_2 = [&]{
		double first{ 0.0 };
		first = 1.0E6*(-D*V_H) / (3.0 * R * T) * K1_(sigma, l[i]);
		first *= integrate(
			[&](double x) { return der_c1(x, expInter[i], slope[i]) * c1(x, expInter[i], slope[i]) / sqrt(x*Pi) / (2.0*x); },
			a[i], b[i]);
		double second = D * integrate(
			[&](double x) { return c1(x, expInter[i], slope[i]) * der2_c1(x, expInter[i], slope[i]); },
			a[i], b[i]);
		double num = first + second;

		double denom{ 0.0 };
		denom = integrate(
			[&](double x) { return sq(c1(x, expInter[i], slope[i]));  },
			a[i], b[i]);
		double val = integrate(
			[&](double x) { return 1.0 / (C * pow((sigma*sqrt(Pi*x)*F_torr(x)), n)); },
			l[i], l[i] + a[i]);
		val /= fff;
		return std::exp(num / denom * val);
	};
	auto criticalConc_formula = [&](int i){
		double tmp{ 0.0 };
		tmp = (K1_(sigma, l[i]) - Kmin) / (Kmax - Kmin);
		tmp = pow(tmp, alpha);
		tmp = pow(1.0 - tmp, 1.0 / beta);
		return tmp / omega;
	};
	double LmaxVal = Lmax();
	cout << LmaxVal << " - max\n";
	while (l[i] < Lmax())
	{
		a[i] = a_func(l[i], delta, a0, l0, mm, Kmax, sigma, alpha2, beta2);
		b[i] = 500.0 * a[i];
		std::tie(expInter[i], slope[i]) = expParams(i);

		m = m_formula();
		B = B_formula();
		deltaT = deltaT_formula();
		data[i][6] = deltaT;

		if (fff != 0.0)
		{
			double Kcur = K1_(sigma, l[i]);
			double sp_fat = Speed_fatigue(Kcur, fff);
			double val = integrate([&](double x) {
				return 1.0 / (C * pow((sigma*sqrt(Pi*x)*F_torr(x)), n));
			}, l[i], l[i] + a[i]);
			val /= fff;

			if (deltaT >= val)
			{
				data[i][6] = val;
				j++;
				if (j == 1)
				{
					data[jjj][4] = i;
					jjj++;
				}
				m = m_formula_2();
				kk = 0;
			}
			else {
				data[i][6] = deltaT;
				kk++;
				if (kk == 1)
				{
					data[jjj][4] = i;
					jjj++;
				}
				j = 0;
			}
		}
		data[i][11] = deltaT;
		data[i][1] = l[i];
		data[1][2] = K1_(sigma, l[i]);
		data[i][3] = m_formula() * omega;
		data[i][4] = deltaT;
		if (fff != 0.0)
		{
			data[i][5] = integrate([&](double x) {
				return 1.0 / (C * pow((sigma*sqrt(Pi*x)*F_torr(x)), n));
			}, l[i], l[i] + a[i]);
			data[i][5] /= fff;
		}
		if (i == 1)
		{
			data[i + 1][7] = data[i][6];
		}
		else if (i > 1)
		{
			data[i + 1][7] = data[i][7] + data[i][6];
		}
		data[i][8] = integrate(
			[&](double x){ return c1(x, expInter[i], slope[i]); }, 0.0, a[i]) / a[i];

		l[i + 1] = l[i] + a[i];
		a[i + 1] = a_func(l[i + 1], delta, a0, l0, mm, Kmax, sigma, alpha2, beta2);
		b[i + 1] = 500.0 * a[i + 1];
		kb[i + 1] = m * c1(a[i] + 500.0*a[i + 1], expInter[i], slope[i]);
		ka[i + 1] = m * c1(a[i + 1] + a[i], expInter[i], slope[i]);
		double newInter, newslope;
		std::tie(newInter, newslope) = expParams(i + 1);
		data[i][9] = a[i];
		data[i][10] = data[i][6];
		data[i][12] = a[i] / data[i][6];
		data[i][13] = expInter[i];
		data[i][14] = slope[i];

		double avConc = integrate(
			[&](double x) { return c1(x, newInter, newslope); }, 0.0, a[i + 1]
			) / a[i + 1];
		double criticalConc = criticalConc_formula(i + 1);
		/*if (K1_(sigma, l[i + 1]) >= Kmax)
			break;*/
		if (deltaT <= 0.0 || avConc >= criticalConc)
			break;
		i++;
	}
	std::cout << i << " iterations\n";
	if (i == 1)
	{
		data[i][14] = 0.0;
		data[1][1] = l0;
		return subMatrix(data, 0, i, 0, 20);
	}
	return subMatrix(data, 0, i - 1, 0, 15);
}

double Fatigue_t_length(double l, double ff, double sigma, double l0)
{
	return integrate([&](double x){ return 1.0 / (C * std::pow(sigma * std::sqrt(Pi*x), n)); }, l0, l, 5.0);
}

std::tuple<Vector, Vector> Fatigue_l_t(double ff, double sigma, double l0, double delta, double Kmax)
{
	auto Lmax = [&] {
		return (1.0 - delta) * sq(Kmax / sigma) / Pi;
	};

	double lMax{ Lmax() };
	int steps{ 100 };
	double step = (lMax - l0) / steps;
	Vector L(steps);

	for (int i{ 0 }; i < steps; ++i)
		L[i] = step * i + l0;

	std::clog << L.size() << ' ' << L.capacity() << '\n';

	Vector T(L.size());
	//std::transform(L.begin(), L.end(), T.begin(), [&](double x) { return Fatigue_t_length(x, ff, sigma, l0); });
	for (int i{ 0 }; i < L.size(); ++i)
		T[i] = Fatigue_t_length(L[i], ff, sigma, l0);


	return std::make_tuple(L, T);
}

double l0_init = 0.003;

vector<std::pair<double, double>> get_Tend_L0_FatMode(double ff, double delta, double Kmax, double sigma, int steps)
{
	using std::pair;
	using std::make_pair;

	auto Lmax = [&] {
		return (1.0 - delta) * sq(Kmax / sigma) / Pi;
	};
	cout << steps << '\n';
	vector<pair<double, double>> plotData;
	double L0_low { l0_init };
	double L0_high { 0.04 };
	double step = (L0_high - L0_low) / double(steps);

	cout << L0_low << ' ' << L0_high << ' ' << step << '\n';

	for (double x{ L0_low }; x <= L0_high; x += step)
	{
		plotData.push_back(make_pair(x, (Fatigue_t_length(Lmax(), ff, sigma, x)) / ff));
	}

	return plotData;
}

void function_calcTime_vs_life()
{

	ofstream outs("out.txt");
	double l0_init = 0.004;
	double l0_step = (0.045 - 0.004) / 50.0;
	for (int i = 0; i < 70; i++)
	{
		if (abs(l0_init + i*l0_step - 0.011) < 0.0001)
			l0_step *= 0.75;
		if (abs(l0_init + i*l0_step - 0.02) < 0.0001)
			l0_step = (0.045 - 0.004) / 50.0;

		auto t0 = high_resolution_clock::now();
		auto res = Solve_env(l0_init + i*l0_step, 0.05, 1.0E-5, 140.0, 10, 80.0, 2, 2, 10.0, 2, 2, 2.5);
		auto t1 = high_resolution_clock::now();
		outs << l0_init + i*l0_step << ' ' << res[res.size() - 1][7] << ' ' << duration_cast<milliseconds>(t1 - t0).count() << '\n';
		cout << i << '\n';
	}

}

std::vector<double> readData(std::istream& is)
{
	std::vector<double> data(12);
	int i{ 0 };

	while (is >> data[i])
		i++;

	if (i != data.size())
	{
		std::cerr << "Wrong input data format: " << data.size() << ' ' << i << '\n';
		exit(4);
	}

	return data;
}

void ProgEnv(int argc, char* argv[])
{
	if (argc < 6) {
		cerr << "Wrong input arguments\n";
		exit(5);
	}

	double l0 = stod(argv[1]);
	double delta = 0.05;
	double a0 = 1.0E-5;

	double sigma = stod(argv[2]);
	double K_min = 10.0;
	double K_max = stod(argv[3]);

	double alpha1 = 2.0;
	double beta1 = 2.0;
	double B = 10.0;

	double alpha2 = 2.0;
	double beta2 = 2.0;
	double Omega = stod(argv[4]);

	auto t0 = high_resolution_clock::now();

	auto res = Solve_env(l0, delta, a0, sigma, K_min, K_max, alpha1, beta1, B, alpha2, beta2, Omega);

	auto t1 = high_resolution_clock::now();

	ofstream outputStream(argv[5]);
	outputStream << duration_cast<milliseconds>(t1 - t0).count() << '\n';

	for (int i = 0; i < res.size(); ++i)
	{
		outputStream << res[i][7] /* t */ << ' ' << res[i][1]  /* l(t) */ << '\n';
	}

}

void Tend_l0(double f)
{
	vector<pair<double, double>> end_plot_data;

	int n = 0;
	for (double l0 = 0.005; l0 <= 0.05; l0 += (0.08 - 0.005) / 40.0)
	{
		ostringstream fileName;
		fileName << n++ << "_.txt";
		ofstream outputStream(fileName.str());

		auto res = Solve_cycle(l0, 0.05, 1.0E-5, 140.0, 10.0, 80.0, 2, 2, 10.0, 2, 2, f, 2.5);

		for (int i = 0; i < res.size(); ++i)
		{
			outputStream << res[i][7] << ' ' << res[i][1] << '\n';
		}
		end_plot_data.push_back(make_pair(l0, res.back()[7]));
	}
}

double life(double l0, double f)
{
	cout << l0 << ' ' << f << '\n';
	//ostringstream fileName;
	//fileName << l0 << '_' << f << "_.txt";
	//ofstream outputStream(fileName.str());

	auto res = Solve_cycle(l0, 0.05, 1.0*1.0E-5, 220.0, 10.0, 80.0, 2, 2, 10.0, 2, 2, f, 1.5);

	//for (int i = 0; i < res.size(); ++i)
	//{
	//	outputStream << res[i][7] << ' ' << res[i][1] << '\n';
	//}
	return res.back()[7];
}

int main(int argc, char* argv[])
{
	//auto res = get_Tend_L0_FatMode(stod(argv[1]), 0.05);

	//return 0;
	cout << argc << '\n';
	l0_init = 0.001;
	if (argc == 2)
	{
		double f = stod(argv[1]);
		ostringstream fileName;

		fileName << f << '_' << "_t_l.txt";

		ofstream endPlotData(fileName.str());
		double step = (0.04 - l0_init) / 30.0;

		for (double l0 = l0_init; l0 <= 0.04; l0 += step)
		{
			endPlotData << l0 << ' ' << life(l0, f) << '\n';
		}
	}
	else {
		double f = stod(argv[1]);
		auto res = get_Tend_L0_FatMode(f, 0.05, 80.0, 220.0, 40);
		ostringstream fileName;
		fileName << f << '_' << "fat_t_l.txt";

		ofstream endPlotData(fileName.str());
		for (const auto p : res)
		{
			endPlotData << p.first << ' ' << p.second << '\n';
		}
	}
}
