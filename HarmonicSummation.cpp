#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\HarmonicSummation.h"
#include <algorithm>

void HarmonicIndex(const std::vector<double>& SegmentFFT, double sampleRate, int nbOfHarmonicsInit, double f0Area, std::vector<int>& harmonics);

double HarmonicSummation(const std::vector<double>& SegmentFFT,  int f0LimitsInf, int f0LimitsSup, int nbOfHarmonicsInit, double sampleRate)
{
	
	std::vector<double> f0Grid;
	std::vector<int>    harmonics(nbOfHarmonicsInit);
	std::vector<double> harmonicsSum;
	double pitch(0);
	double sumHarmo(0);
	double maxSumHarmo(0);

	for (double i = f0LimitsInf; i < f0LimitsSup; i += 0.1)
	{
		HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, i, harmonics);
		sumHarmo = 0;

		for (int j = 0; j < harmonics.size(); ++j)
		{
			sumHarmo += SegmentFFT[(harmonics[j])];
		}

		if (maxSumHarmo < sumHarmo)
		{
			pitch = i;
			maxSumHarmo = sumHarmo;
		}
	}
	
	//for each frequence we get the sum of the value of the exact harmonic.
	//make it more precise :
	/*
	for (double i = pitch-4.0; i < pitch+4.0; i += 0.01)
	{
		HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, i, harmonics);
		sumHarmo = 0.0;

		for (int j = 0; j < harmonics.size(); ++j)
		{
			sumHarmo += SegmentFFT[(harmonics[j])];
		}

		if (maxSumHarmo < sumHarmo)
		{
			pitch = i;
			maxSumHarmo = sumHarmo;
		}
	}

	for (double i = pitch - 0.1; i < pitch + 0.1; i += 0.001)
	{
		HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, i, harmonics);
		sumHarmo = 0.0;

		for (int j = 0; j < harmonics.size(); ++j)
		{
			sumHarmo += SegmentFFT[(harmonics[j])];
		}

		if (maxSumHarmo < sumHarmo)
		{
			pitch = i;
			maxSumHarmo = sumHarmo;
		}
	}*/
	return(pitch);
}

void HarmonicIndex(const std::vector<double>&  SegmentFFT, double sampleRate, int nbOfHarmonicsInit, double f0Area, std::vector<int>& harmonics)
{
	size_t NFFT = SegmentFFT.size();

	for (int i = 1; i < nbOfHarmonicsInit+1; ++i)
	{
		harmonics[(__int64)i - (__int64)1] = static_cast<int>(round((i * f0Area) * (2.0 * NFFT / sampleRate)));//*2 to get 2^19
	}
	
}