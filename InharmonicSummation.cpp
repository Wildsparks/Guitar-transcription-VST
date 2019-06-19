#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\InharmonicSummation.h"
#include <algorithm>
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Source\package eigen\Eigen\Dense"
#include <iostream>

std::vector<int> inharmonicIndex(std::vector<float>& const SegmentFFT, unsigned int sampleRate, unsigned int nBOfHarmonic, double pitch, double beta);

void InharmonicSummation(std::vector<float>& const SegmentFFT, float pitchInitial, unsigned int highNbOfHarmonics, unsigned int sampleRate, double maxBetaGrid, double minBetaGrid, float betaRes, float lengthFFT, float& pitchEstimate, float& BEstimate, float& costFunctionMaxVal)
{
	std::vector<float> upperLowerOctave(2);
	std::vector<double> maxCost(2);
	std::vector<double> pitchCandidate(2);
	std::vector<double> betaCandidate(2);

	upperLowerOctave[1] = (2 * pitchInitial);
	upperLowerOctave[0] = (pitchInitial);

	std::vector<double> pitchGrid;
	double pitchWidth ( (3.0 + 3.0) * lengthFFT / pow(2,18) * sampleRate / lengthFFT);
	int counterBeta(0);

	int betaGridSize(ceil(double((maxBetaGrid - minBetaGrid) / betaRes)));

	for(int octave=0; octave<2; ++octave)
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
		betaCandidate[octave] = double(indexBeta) * betaRes+ minBetaGrid;

	}//done for lower and upper octave
	
	int maxIndexOctave(0);
	int maxvalueOctave(0);

	for (int i = 0; i < 2; ++i)
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
std::vector<int> inharmonicIndex(std::vector<float>& const SegmentFFT, unsigned int sampleRate, unsigned int nBOfHarmonic, double pitch, double beta)
{
	std::vector<double> phi(nBOfHarmonic);
	std::vector<int> index(nBOfHarmonic);

	phi[0] = pitch;
	index[0] = round(phi[0] * (2 * SegmentFFT.size() / sampleRate) + 1);

	for (int j = 1; j < nBOfHarmonic; ++j)
	{
		phi[j] = pitch * j * sqrt(1 + beta * pow(j, 2));
		index[j] = round(phi[j] * (2 * SegmentFFT.size() / sampleRate) + 1);
	}

	return(index);
}
