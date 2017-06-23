#ifndef CRACK_MODEL_H
#define CRACK_MODEL_H

#include <cmath>
#include "Support.hpp"


class CrackGrowthModel
{
public:
	double l;
	double delta;
	double a0;
	double l0;
	double mm;
	double K1c_max;
	double sigma;
	double alpha_2, beta_2;

	void nextState()
	{
		state.l += prerupAreaLen(state);
	}

	double prerupAreaLen()
	{
		double L_max = L_max();
		double Alpha = std::pow((L_max - l) / (L_max - l0), beta_2);
		Alpha = std::pow(1.0 - Alpha, alpha_2);

		return a0*(Alpha*(B - 1.0) + 1.0);
	}

	double L_max()
	{
		return sq(K_max) / (Pi * sigma)
	}
};

#endif // CRACK_MODEL_H