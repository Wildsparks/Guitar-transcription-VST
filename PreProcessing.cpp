#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"

#include <iostream>

#include <complex>
#include <valarray>
#include <algorithm> 

#include <vector>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

#define M_PI 3.141592653589793238460

void fftprec(CArray& x);
void Hanningprec(int lenth, std::vector<float>& hann);
void ifft(CArray& x);

void Preprocessing(std::vector<float>& segment, int sampleRate)
{

	//=====================================================//
	//==================Hilbert transform==================//
	//=====================================================//
	if (segment.size() % 2 == 1)
	{
		segment.pop_back();
	}

	CArray segmentFreq(segment.size());
	//windows//
	std::vector<float> hann(segment.size());
	Hanningprec(segment.size(), hann);
	//end - windows//

	for (int i = 0;i < segment.size();++i)
	{
		segmentFreq[i] = segment[i]* hann[i];
	}

	fftprec(segmentFreq);

	std::vector<int> h(segment.size());

	segmentFreq[0] = 1.0 * segmentFreq[0];
	for (int i = 1; i < segment.size()/2; ++i)
	{
		segmentFreq[i] = 2.0 * segmentFreq[i];
	}
	for (int i = segment.size() / 2; i < segment.size(); ++i)
	{
		segmentFreq[i] = 1.0 * segmentFreq[i];
	}

	ifft(segmentFreq);

	for (int i = 0;i < segment.size();++i)
	{
		segment[i] = segmentFreq[i].real();
	}


	//=====================================================//
	//===============Hilbert transform end=================//
	//=====================================================//
}













void fftprec(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fftprec(even);
	fftprec(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

void Hanningprec(int lenth, std::vector<float>& hann)
{
	//Hanning Window
	for (int i = 0; i < lenth; i++) {
		double multiplier = 0.5 * (1 - cos(2 * M_PI * i / (lenth - 1.0)));
		hann[i] = multiplier;
	}
	//
}

void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fftprec(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}