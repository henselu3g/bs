// BlackSholes.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <string>
//#include <random>
#include <math.h>

double normal_s_CFD(double value)
{
	return 0.5 * erfc(-value * M_SQRT1_2);
}

double normal_s_PDF(double value)
{
	return (.5 * M_2_SQRTPI * M_SQRT1_2) * exp(-0.5 * value * value);
}

class BlackScholes
{
public:
//private:
	double S;
	double K;
	double rf;
	double div_rate;
	double vol;
	double yrs_to_maturity;

	double c_pr;
	double c_delta;
	double c_gamma;
	double c_theta;
	double c_theta_per_day;

//public:
	BlackScholes() { };
	BlackScholes(double s, double k, double r, double d, double v, double yrs);

	~BlackScholes() {};
};

BlackScholes::BlackScholes(double s, double k, double r, double d, double v, double yrs) :
	S(s), K(k), rf(r), div_rate(d), vol(v), yrs_to_maturity(yrs)
{ 
	double d1 = (log(S / K) + (rf - div_rate + .5 * vol * vol) * yrs_to_maturity) / (vol * sqrt(yrs_to_maturity));
	double d2 = (log(S / K) + (rf - div_rate - .5 * vol * vol) * yrs_to_maturity) / (vol * sqrt(yrs_to_maturity));

	double nd1 = normal_s_CFD(d1);
	double nd2 = normal_s_CFD(d2);

	double nd1p = normal_s_PDF(d1);
	double nd2p = normal_s_PDF(d2);

	c_pr = S*exp(-div_rate * yrs_to_maturity)*nd1 - K * exp(-rf*yrs_to_maturity)*nd2;
	c_delta = exp(-div_rate * yrs_to_maturity) * nd1;
	c_gamma = exp(-div_rate * yrs_to_maturity) * nd1p / (S *vol * sqrt(yrs_to_maturity));
	c_theta = -0.5 * vol * S * exp(-div_rate * yrs_to_maturity) * nd1p / sqrt(yrs_to_maturity) + div_rate * S * exp(-div_rate * yrs_to_maturity) * nd1 - rf * K * exp(-rf*yrs_to_maturity)*nd2;
	c_theta_per_day = c_theta / 365;
};

int main()
{
	double S = 40;
	double K = 40;
	double vol = .3;
	double rf = .08;
	double div_rate = 0;
	int days_to_maturity = 91;
	double yrs_to_maturity = static_cast<double>(days_to_maturity) / 365;

	BlackScholes bs0(S, K, rf, div_rate, vol, yrs_to_maturity);
	double call_pr_40 = bs0.c_pr;

	S = 20;
	while (S < 60.5) {
		BlackScholes bs(S, K, rf, div_rate, vol, yrs_to_maturity);
		BlackScholes bs_T(S, K, rf, div_rate, vol, 1./365);
		double profit = call_pr_40 - bs.c_pr;
		double profit_T = call_pr_40 - bs_T.c_pr;



		S += 0.5;
	}
	
	return 0;
}

