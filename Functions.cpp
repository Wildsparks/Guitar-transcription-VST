#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Functions.h"

void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

void Hanning(int lenth, std::vector<double>& hann)
{
	//Hanning Window
	for (int i = 0; i < lenth; i++) {
		double multiplier = 0.5 * (1.0 - cos(2.0 * M_PI * i / (lenth - 1.0)));
		hann[i] = multiplier;
	}
	//
}

void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= (double)x.size();
}

void Gaussian(int lenth, std::vector<double>& Gaussian)
{
	double Nwin = lenth-1.0;
	double a = 2.5;

	for (double i = -Nwin /2.0; i < Nwin /2.0; ++i)
	{
		Gaussian[unsigned int(double(i) + double(Nwin / 2.0))]= exp(-(1.0 / 2.0) * pow((a * i / (Nwin / 2.0)),2.0));
	}
}
