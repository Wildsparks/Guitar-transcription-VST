#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\AmplitudesEstimation.h"

#include <complex>
#include <cmath>
#include <iostream>


#define M_PI 3.14159265359

using namespace std;

std::vector<double> AmplitudesEstimation(std::vector<std::complex<double>> const& segment, double f0Features, double betaFeatures, int lengthOfSegment, double sampleRate, int highNbOfHarmonics)
{

	Eigen::MatrixXcf Z(lengthOfSegment, highNbOfHarmonics);
	const complex<double> i(0.0, 1.0);

	for (int m = 1; m < highNbOfHarmonics+1; ++m)
	{
		for (int n = 0; n < lengthOfSegment; ++n)
		{
			Z(n, m-1) = std::exp((i * 2.0 * M_PI * f0Features * double(m) * std::sqrt(1.0 + betaFeatures * pow(m, 2.0))) / sampleRate * double(n));
		}
	}

	Eigen::VectorXcf alpha(highNbOfHarmonics);
	Eigen::VectorXcf x(lengthOfSegment);
	
	for (int c = 0; c < lengthOfSegment; ++c)
	{
		x(c) = segment[c];
	}

	Eigen::MatrixXcf ZT((Z.transpose()).conjugate());
	Eigen::MatrixXcf ZTINV((ZT * Z).inverse());

	alpha = ZTINV * ZT * x;

	std::vector<double> amplitudesAbs(highNbOfHarmonics);
	for (int c = 0; c < highNbOfHarmonics; ++c)
	{
		amplitudesAbs[c] = std::abs(alpha(c));
	}

	return(amplitudesAbs);
}
