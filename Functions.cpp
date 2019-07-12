#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Functions.h"

void Hanning(int lenth, std::vector<double>& hann)
{
	//Hanning Window
	for (int i = 0; i < lenth; i++) {
		double multiplier = 0.5 * (1.0 - cos(2.0 * M_PI * i / (lenth - 1.0)));
		hann[i] = multiplier;
	}
	//
}

void Gaussian(int lenth, std::vector<double>& Gaussian)
{
	double Nwin = lenth-1.0;
	double a = 2.5;

	for (double i = -Nwin /2.0; i < Nwin /2.0; ++i)
	{
		Gaussian[static_cast<__int64>(i + (Nwin / 2.0))]= exp(-(1.0 / 2.0) * pow((a * i / (Nwin / 2.0)),2.0));
	}
}
