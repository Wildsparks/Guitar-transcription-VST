
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PluckingPositionEstimatorLSD.h"

#include <algorithm>

#define M_PI 3.141592653589793238460

float PluckingPositionEstimatorLSD(std::vector<float> amplitudesAbs, double L)
{
	double delta ( 1 );
	double M = amplitudesAbs.size();

	std::vector<double> Cm(M);
	std::vector<double> D_LS(0);

	double P(0);
	double sumAmpAbsOverCm;

	D_LS.clear();

	for (double cnt = 6; cnt < ceil(L / 2); cnt += 0.001)
	{
		
		P = cnt / L;

		for (int m = 0; m < M; ++m)
		{
			Cm[m] = 2 * delta / (pow(m+1.0 , 2) * pow(M_PI , 2) * P * (1 - P))*abs(sin(((m+1.0)*M_PI * P)));
		}

		sumAmpAbsOverCm = 0;

		for (int m = 0; m < M; ++m)
		{
			sumAmpAbsOverCm += pow((10 * log10(amplitudesAbs[m] / Cm[m])), 2);
		}

		D_LS.push_back( sqrt(1.0 / M * sumAmpAbsOverCm) );

	}

	int minElementIndex = std::min_element(D_LS.begin(), D_LS.end()) - D_LS.begin();

	return(double(minElementIndex)*0.001+6.0);

}