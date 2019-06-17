#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\HarmonicSummation.h"
#include <algorithm>
#define f0GridStep 0.1

std::vector<float> HarmonicIndex(std::vector<float> SegmentFFT, unsigned int sampleRate, unsigned int nbOfHarmonicsInit, double f0Area);

float HarmonicSummation(std::vector<float> SegmentFFT,  unsigned int f0LimitsInf, unsigned int f0LimitsSup, unsigned int nbOfHarmonicsInit, unsigned int sampleRate)
{
	
	std::vector<float> f0Grid;
	std::vector<float> harmonics;
	std::vector<float> harmonicsSum;
	float pitch;

	for (double i = f0LimitsInf; i < f0LimitsSup; i += 0.1)
	{
		f0Grid.push_back(i);
	}

	harmonicsSum.resize(f0Grid.size(),0);

	for (int i = 0; i < f0Grid.size(); i++)
	{
		harmonics=HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, f0Grid[i]);

		for (int j = 0; j < harmonics.size(); ++j)
		{
			harmonicsSum[i] += SegmentFFT[int(harmonics[j])];
		}
		
	}

	pitch = f0Grid[std::max_element(harmonicsSum.begin(), harmonicsSum.end()) - harmonicsSum.begin()];

	/*//for each frequence we get the sum of the value of the exact harmonic.
	//make it more precise :
	f0Grid.clear();
	harmonics.clear();
	harmonicsSum.clear();

	for (double i = pitch-4.0; i < pitch+4.0; i += 0.01)
	{
		f0Grid.push_back(i);
	}

	harmonicsSum.resize(f0Grid.size(), 0);

	for (int i = 0; i < f0Grid.size(); i++)
	{
		harmonics = HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, f0Grid[i]);

		for (int j = 0; j < harmonics.size(); ++j)
		{
			harmonicsSum[i] += SegmentFFT[int(harmonics[j])];
		}

	}

	pitch = f0Grid[std::max_element(harmonicsSum.begin(), harmonicsSum.end()) - harmonicsSum.begin()];

	//for each frequence we get the sum of the value of the exact harmonic.
	//make it more precise :
	f0Grid.clear();
	harmonics.clear();
	harmonicsSum.clear();

	for (double i = pitch - 0.01; i < pitch + 0.01; i += 0.001)
	{
		f0Grid.push_back(i);
	}

	harmonicsSum.resize(f0Grid.size(), 0);

	for (int i = 0; i < f0Grid.size(); i++)
	{
		harmonics = HarmonicIndex(SegmentFFT, sampleRate, nbOfHarmonicsInit, f0Grid[i]);

		for (int j = 0; j < harmonics.size(); ++j)
		{
			harmonicsSum[i] += SegmentFFT[int(harmonics[j])];
		}

	}

	pitch = f0Grid[std::max_element(harmonicsSum.begin(), harmonicsSum.end()) - harmonicsSum.begin()];*/

	/*//return maximum near f0 estimate
	std::vector<float> lowerIndex= HarmonicIndex(SegmentFFT, sampleRate, 1, pitch - 0.75 * pitch);
	std::vector<float> upperIndex= HarmonicIndex(SegmentFFT, sampleRate, 1, pitch + 0.75 * pitch);

	pitch = *std::max_element(SegmentFFT.begin() + lowerIndex[0], SegmentFFT.begin() + upperIndex[0])*/

	return(pitch);

}


std::vector<float> HarmonicIndex(std::vector<float> SegmentFFT, unsigned int sampleRate, unsigned int nbOfHarmonicsInit, double f0Area)
{
	std::vector<float> harmonics;
	double NFFT = SegmentFFT.size();

	for (int i = 1; i < nbOfHarmonicsInit; ++i)
	{
		harmonics.push_back(round((i*f0Area)*(2*NFFT/sampleRate)));
	}
	
	return(harmonics);
}