
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PluckingPositionEstimatorLSD.h"

#include <algorithm>

#define M_PI 3.141592653589793238460

double PluckingPositionEstimatorLSD(std::vector<double> amplitudesAbs, double L)
{
	double delta ( 1.0 );
	double M = amplitudesAbs.size();

	std::vector<double> Cm(M);
	//std::vector<double> D_LS(0);

	double P(0);
	double sumAmpAbsOverCm;
	double minValue = std::numeric_limits<double>::max();
	double result(0);

	//D_LS.clear();

	for (double cnt = 6.0; cnt < ceil(L / 2.0); cnt += 0.001) // search grid in cm from the bridge.
	{
		
		P = cnt / L;

		for (int m = 0; m < M; ++m)
		{
			Cm[m] = 2.0 * delta / (pow(m+1.0 , 2.0) * pow(M_PI , 2.0) * P * (1.0 - P))*std::abs(sin(((m+1.0)*M_PI * P)));
		}

		sumAmpAbsOverCm = 0;

		for (int m = 0; m < M; ++m)
		{
			sumAmpAbsOverCm += pow((10.0 * log10(amplitudesAbs[m] / Cm[m])), 2.0);
		}

		//D_LS.push_back( sqrt(1.0 / M * sumAmpAbsOverCm) );
		if ((sqrt(1.0 / M * sumAmpAbsOverCm)) < minValue)
		{
			minValue = sqrt(1.0 / M * sumAmpAbsOverCm);
			result = cnt;
		}
	}

	//int minElementIndex = std::min_element(D_LS.begin(), D_LS.end()) - D_LS.begin();
	return(result);
	//return(double(minElementIndex)*0.001+6.0);

}