#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\InharmonicSummation.h"
#include <algorithm>
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Source\package eigen\Eigen\Dense"
#include <iostream>

std::vector<double> inharmonicIndex(const std::vector<double>& SegmentFFT, double sampleRate, int nBOfHarmonic, double pitch, double beta);

void InharmonicSummation(const std::vector<double>& SegmentFFT, double pitchInitial, int highNbOfHarmonics, double sampleRate, double maxBetaGrid, double minBetaGrid, double betaRes, double lengthFFT, double& pitchEstimate, double& BEstimate, double& costFunctionMaxVal)
{
	std::vector<double> upperLowerOctave(2);
	std::vector<double> maxCost(2);
	std::vector<double> pitchCandidate(2);
	std::vector<double> betaCandidate(2);

	upperLowerOctave[1] = (2 * pitchInitial);
	upperLowerOctave[0] = (pitchInitial);

	std::vector<double> pitchGrid;
	double pitchWidth ( (3.0+3.0) / pow(2.0,18.0) * sampleRate );

	int counterBeta(0);

	int betaGridSize(static_cast<int>(ceil(((maxBetaGrid) - (minBetaGrid)) / betaRes)));

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
				std::vector < double > index = inharmonicIndex(SegmentFFT, sampleRate, highNbOfHarmonics, pitchGrid[k], i);

				while (*(index.end()-1) > SegmentFFT.size())
				{
					index.pop_back();
				}

				costFunction(counterBeta, k) = 0;

				for (int m = 0; m < index.size(); ++m)
				{
					if(index[m]>=0)
					{
						costFunction(counterBeta, k) += SegmentFFT[index[m]] * pow(m, 2.0);
					}
					else
					{
						costFunction(counterBeta, k) += 0;
					}
					
				}
			}
			counterBeta++;
			
		}

		maxCost[octave]=(costFunction(0, 0));
		int indexBeta(0);
		int indexPitch(0);

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
	double maxvalueOctave(-1.0);

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
std::vector<double> inharmonicIndex(const std::vector<double>& SegmentFFT, double sampleRate, int nBOfHarmonic, double pitch, double beta)
{
	std::vector<double> phi  (nBOfHarmonic + double(1));
	std::vector<double> index(nBOfHarmonic + double(1));

	phi  [0] = pitch;
	index[0] = static_cast<int>(round(phi[0] * (2.0 * SegmentFFT.size() / sampleRate)));

	for (int j = 1; j < (nBOfHarmonic+1); ++j)
	{
		phi[j]   = pitch * (j+1.0) * sqrt(1.0 + beta * pow((j+1.0), 2.0));
		index[j] = (round(phi[j] * (2.0 * SegmentFFT.size() / sampleRate)));
	}

	return(index);
}
