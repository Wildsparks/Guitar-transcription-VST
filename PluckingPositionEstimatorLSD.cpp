
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PluckingPositionEstimatorLSD.h"

#include <algorithm>

#define M_PI 3.141592653589793238460

double PluckingPositionEstimatorLSD(std::vector<double> amplitudesAbs, double L)
{
	double delta ( 1.0 );
	double M = amplitudesAbs.size();

	std::vector<double> Cm(M);
	std::vector<double> diffvector(M-1);

	double P(0);
	double J;
	double minValue = std::numeric_limits<double>::max();
	double result(0);

	for (double cnt = 6.0; cnt < ceil(L / 2.0); cnt += 0.1) // search grid in cm from the bridge.
	{
		
		P = cnt / L;

		for (int m = 0; m < M; ++m)
		{
			Cm[m] = 2.0 * delta / (pow(m+1.0 , 2.0) * pow(M_PI , 2.0) * P * (1.0 - P))*std::abs(sin(((m+1.0)*M_PI * P)));
		}

		for (int m = 1; m < M; ++m)
		{
			diffvector[m-1] = (log10( std::abs( Cm[m] ) ) - log10(std::abs(Cm[m-1]))) - (log10(std::abs(amplitudesAbs[m])) - log10(std::abs(amplitudesAbs[m - 1])));
		}

		J = 0;

		for (int m = 0; m < M-1; ++m)
		{
			J += pow((std::abs(diffvector[m])), 2.0);
		}

		if (J < minValue)
		{
			minValue = J;
			result = cnt;
		}
	}

	return(result);

}

double PluckingPositionEstimatorLSDold(std::vector<double> amplitudesAbs, double L)
{
	double delta(1.0);
	double M = amplitudesAbs.size();

	std::vector<double> Cm(M);

	double P(0);
	double sumAmpAbsOverCm;
	double minValue = std::numeric_limits<double>::max();
	double result(0);

	for (double cnt = 6.0; cnt < ceil(L / 2.0); cnt += 0.1) // search grid in cm from the bridge.
	{

		P = cnt / L;

		for (int m = 0; m < M; ++m)
		{
			Cm[m] = 2.0 * delta / (pow(m + 1.0, 2.0) * pow(M_PI, 2.0) * P * (1.0 - P)) * std::abs(sin(((m + 1.0) * M_PI * P)));
		}

		sumAmpAbsOverCm = 0;

		for (int m = 0; m < M; ++m)
		{
			sumAmpAbsOverCm += pow((10.0 * log10(amplitudesAbs[m] / Cm[m])), 2.0);
		}

		if ((sqrt(1.0 / M * sumAmpAbsOverCm)) < minValue)
		{
			minValue = sqrt(1.0 / M * sumAmpAbsOverCm);
			result = cnt;
		}
	}

	return(result);

}