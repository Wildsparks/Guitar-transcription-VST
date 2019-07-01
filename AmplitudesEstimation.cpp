#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\AmplitudesEstimation.h"

#include <complex>
#include <cmath>

#define M_PI 3.141592653589793238460

using namespace std;

std::vector<float> AmplitudesEstimation(std::vector<float> const& segment, double f0Features, double betaFeatures, int lengthOfSegment, unsigned int sampleRate, unsigned int highNbOfHarmonics)
{
	Eigen::MatrixXcd Z(lengthOfSegment, highNbOfHarmonics);
	const   complex<double> i(0.0, 1.0);
	complex<double> plop;

	for (double m = 0; m < highNbOfHarmonics; ++m)
	{
		for (double n = 0; n < lengthOfSegment; ++n)
		{
			Z(n, m) = exp((i * 2.0 * M_PI * f0Features * m * sqrt(1.0 + betaFeatures * pow(m, 2.0))) / double(sampleRate) * (n)); //I have delete a "-1"
		}
	}

	Eigen::VectorXcd alpha(highNbOfHarmonics);
	Eigen::VectorXd x(lengthOfSegment);
	std::vector<float> amplitudesAbs(highNbOfHarmonics);
	
	for (int c = 0; c < lengthOfSegment; ++c)
	{
		x(c) = segment[c];
	}
	
	alpha = ((Z.transpose() * Z).inverse() * Z.transpose()) * x;
	
	for (int c = 0; c < highNbOfHarmonics; ++c)
	{
		amplitudesAbs[c] = abs(alpha(c));
	}

	return(amplitudesAbs);
}
