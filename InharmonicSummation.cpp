#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\InharmonicSummation.h"
#include <algorithm>
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Source\package eigen\Eigen\Dense"
#include <iostream>

std::vector<int> inharmonicIndex(const std::vector<float>& SegmentFFT, unsigned int sampleRate, unsigned int nBOfHarmonic, double pitch, double beta);

void InharmonicSummation(const std::vector<float>& SegmentFFT, float pitchInitial, unsigned int highNbOfHarmonics, unsigned int sampleRate, double maxBetaGrid, double minBetaGrid, double betaRes, double lengthFFT, float& pitchEstimate, float& BEstimate, float& costFunctionMaxVal)
{
	std::vector<float> upperLowerOctave(2);
	std::vector<double> maxCost(2);
	std::vector<double> pitchCandidate(2);
	std::vector<double> betaCandidate(2);

	upperLowerOctave[1] = (2 * pitchInitial);
	upperLowerOctave[0] = (pitchInitial);

	std::vector<double> pitchGrid;
	double pitchWidth ( (3.0 + 3.0) / pow(2.0,18.0) * sampleRate );

	int counterBeta(0);

	int betaGridSize(ceil(double((maxBetaGrid - minBetaGrid) / betaRes)));

	for(int octave=0; octave<1; ++octave)//2 if upper octave
	{ 
		for (double j = upperLowerOctave[octave] - pitchWidth; j < upperLowerOctave[octave] + pitchWidth; j += sampleRate / lengthFFT)
		{
			pitchGrid.push_back(j);
		}

		Eigen::MatrixXd costFunction(betaGridSize, pitchGrid.size());
		counterBeta = 0;

		for (double i = minBetaGrid; i < maxBetaGrid; i += betaRes)
		{
			for (int k = 0; k < pitchGrid.size(); ++k)
			{
				std::vector<int> index = inharmonicIndex(SegmentFFT, sampleRate, highNbOfHarmonics, pitchGrid[k], i);

				while (*(index.end()-1) > SegmentFFT.size())
				{
					index.pop_back();
				}

				costFunction(counterBeta, k) = 0;

				for (int m = 0; m < index.size(); ++m)
				{
					costFunction(counterBeta, k) += SegmentFFT[index[m]];
				}
			}
			counterBeta++;
			
		}

		maxCost[octave]=(costFunction(0, 0));
		unsigned int indexBeta(0);
		unsigned int indexPitch(0);

		for (int i = 0; i < betaGridSize; ++i)
		{
			for (int j = 0; j < pitchGrid.size(); ++j)
			{
				if (maxCost[octave] < costFunction(i, j))
				{
					maxCost[octave] = costFunction(i, j);
					indexBeta = i;
					indexPitch = j;
				}
			}
		}

		//pitch//
		pitchCandidate[octave] = pitchGrid[indexPitch];
		//beta//
		betaCandidate[octave] = double(indexBeta) * betaRes + minBetaGrid;

	}//done for lower and upper octave
	
	int maxIndexOctave(0);
	float maxvalueOctave(-1.0);

	for (int i = 0; i < 1; ++i)//2 if upper octave
	{
		if (maxCost[i] > maxvalueOctave)
		{
			maxvalueOctave = maxCost[i];
			maxIndexOctave = i;
		}
	}

	pitchEstimate      = pitchCandidate[maxIndexOctave];
	BEstimate          = betaCandidate[maxIndexOctave];
	costFunctionMaxVal = maxvalueOctave;

}
/*
% This function creates an index vector for signal partials
% used for harmonic summation.
% The fundamental frequency input vectorand the FFT - length,
% will create the dimensions of the output.
% The peaks are placed in a zero vector, on the correct frequency axis.
*/
std::vector<int> inharmonicIndex(const std::vector<float>& SegmentFFT, unsigned int sampleRate, unsigned int nBOfHarmonic, double pitch, double beta)
{
	std::vector<double> phi((long)nBOfHarmonic + 1.0);
	std::vector<int> index((long)nBOfHarmonic + 1.0);

	phi[0] = pitch;
	index[0] = round(phi[0] * (2.0 * SegmentFFT.size() / sampleRate));

	for (unsigned int j = 1; j < (nBOfHarmonic+1); ++j)
	{
		phi[j] = pitch * (j+1.0) * sqrt(1.0 + beta * pow((j+1.0), 2.0));
		index[j] = round(phi[j] * (2.0 * SegmentFFT.size() / sampleRate));
	}

	return(index);
}
