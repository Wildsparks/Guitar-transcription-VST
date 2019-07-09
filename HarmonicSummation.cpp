#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\HarmonicSummation.h"
#include <algorithm>

void HarmonicIndex(const std::vector<float>& SegmentFFT, double sampleRate, int nbOfHarmonicsInit, double f0Area, std::vector<int>& harmonics);

double HarmonicSummation(const std::vector<float>& SegmentFFT,  int f0LimitsInf, int f0LimitsSup, int nbOfHarmonicsInit, double sampleRate)
{
	
	std::vector<float> f0Grid;
	std::vector<int>   harmonics(nbOfHarmonicsInit);
	std::vector<float> harmonicsSum;
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
	}

	/*
	//for each frequence we get the sum of the value of the exact harmonic.
	//make it more precise :
	f0Grid.clear();
	harmonicsSum.clear();

	for (double i = pitch - 0.01; i < pitch + 0.01; i += 0.001)
	{
		f0Grid.push_back(i);
	}

	harmonicsSum.resize(f0Grid.size(), 0);

	for (int i = 0; i < f0Grid.size(); i++)
	{
		HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, f0Grid[i], harmonics);

		for (int j = 0; j < harmonics.size(); ++j)
		{
			harmonicsSum[i] += SegmentFFT[(harmonics[j])];
		}

	}

	pitch = f0Grid[std::max_element(harmonicsSum.begin(), harmonicsSum.end()) - harmonicsSum.begin()-1.0];
	/*
	//return maximum near f0 estimate
	std::vector<int> lowerIndex= HarmonicIndex(SegmentFFT, sampleRate, 1, pitch - 0.75 * pitch, harmonics);
	std::vector<int> upperIndex= HarmonicIndex(SegmentFFT, sampleRate, 1, pitch + 0.75 * pitch, harmonics);
	pitch = *std::max_element(SegmentFFT.begin() + lowerIndex[0], SegmentFFT.begin() + upperIndex[0])
	*/
	return(pitch);
}

void HarmonicIndex(const std::vector<float>&  SegmentFFT, double sampleRate, int nbOfHarmonicsInit, double f0Area, std::vector<int>& harmonics)
{
	size_t NFFT = SegmentFFT.size();

	for (unsigned int i = 1; i < nbOfHarmonicsInit+1; ++i)
	{
		harmonics[i - 1] = (round((i * f0Area) * (2.0 * NFFT / sampleRate)));//*2 to get 2^19
	}
	
}