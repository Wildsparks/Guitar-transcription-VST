#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"
#include <iostream>

#include <complex>
#include <valarray>
#include <algorithm> 

#include <vector>
#include <numeric>
#include <string>
#include <functional>

#define TIMEBUFFER 0.04 //40ms of observation
#define M_PI 3.141592653589793238460

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void onsetDetector(std::vector<float> const& storage, std::vector<double>& timeOfOnset, int sampleRate, int const& thresholdValue);
void fft(CArray& x);
void Hanning(int lenth, std::vector<float>& hann);

void Onset(const float* channelData, int bufsize,int sampleRate, std::vector<double>& timeOfOnset, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast, int const& thresholdValue)
{
	float nbEchentillon ((TIMEBUFFER * sampleRate));
	int nbBuffers    (ceil(nbEchentillon/ bufsize));                // number of buffer to store to get TIMEBUFFER miliseconde (ceil get the int sup)
	
	//set the first vector of storage//
	if (storagePast.size() != int(double(bufsize) * nbBuffers * 3)) //120 milisecond
	{

		storagePast.clear();
		storageActual.clear();
		
		while (storagePast.size() != int(double(bufsize) * nbBuffers * 3))
		{
			storagePast.push_back(1);
		}

	}
	//done//

	if (counter == nbBuffers)
	{

		counter = 0;

		for (int i = 0; i < (bufsize * nbBuffers); i++)
		{
			storagePast.push_back(storageActual[i]); //add the new 40ms to the old 80ms 
		}

		//120ms analysis from 20ms to 60 ms to get 40ms after and pluck excitation
		onsetDetector(storagePast, timeOfOnset, sampleRate, thresholdValue);

		for (int i = 0; i < (bufsize * nbBuffers * 3); i++)
		{
			storagePast[i] = (storagePast[i + double(bufsize) * nbBuffers]); //40ms first milisecond to trash by decay
		}

		for (int i = 0; i < (bufsize * nbBuffers); i++)
		{
			storagePast.pop_back(); //delete last 40ms which are copy
		}

		storageActual.clear();
		
	}
	else
	{
		for (int i = 0; i < bufsize; i++)
		{
			storageActual.push_back(channelData[i]); //40ms
		}
		counter++;
	}
}

void onsetDetector(std::vector<float> const& storage, std::vector<double>& timeOfOnset, int sampleRate, int const& thresholdValue)
{

	unsigned int gapBetweenFrame(1);  // in ms
	unsigned int lenthOfFrame(40); // in ms
	unsigned int gapInSample(gapBetweenFrame * pow(10, -3) * double(sampleRate));
	unsigned int nbFrame((float(storage.size()) / float(gapInSample)));
	unsigned int lenthFft(lenthOfFrame * pow(10, -3) * sampleRate);

	unsigned int beginFrame(0 * pow(10, -3) * sampleRate / gapInSample);
	unsigned int stopFrame((0 * pow(10, -3) * sampleRate / gapInSample) + lenthFft / gapInSample);

	//windows//
	std::vector<float> hann(lenthFft);
	Hanning(lenthFft, hann);
	//end - windows//

	CArray frame(lenthFft);
	std::vector <float> frameAbs(lenthFft);
	std::vector <std::vector <float>> spectrogram;

	//Spectrogram//
	for (int i = beginFrame; i < (nbFrame - stopFrame); i++)
	{
		for (int j = 0; j < lenthFft; j++)
		{
			frame[j] = (storage[(i * floor(gapInSample)) + j] * double(hann[j]));
		}

		fft(frame);

		for (int k = 0; k < (lenthFft); k++)
		{
			frameAbs[k] = std::abs(frame[k]);
		}

		spectrogram.push_back(frameAbs);
	}

	//end - Spectrogram//

	float gamma(0.94); // as in the paper in %
	float J(gamma * double(lenthFft));
	std::vector <float> invSparsity;

	float sumFreq2;
	float sumFreq4;

	//inverse sparsity//
	for (int i = 0; i < (nbFrame - stopFrame - beginFrame); i++)
	{
		sumFreq2 = 0;
		sumFreq4 = 0;
		for (int j = 0; j < J; j++)
		{
			spectrogram[i][j] = pow(spectrogram[i][j], 2);
			sumFreq2 += spectrogram[i][j]; //sum(x^2)
			spectrogram[i][j] = pow(spectrogram[i][j], 2);
			sumFreq4 += spectrogram[i][j]; //sum(x^4)
		}
		if (sumFreq4 == 0)
		{
			invSparsity.push_back(0);
		}
		else
		{
			invSparsity.push_back(sumFreq2 / (pow(sumFreq4, 0.25) * pow(J, 0.25)));
		}

	}

	//filter if it's necessary
	/*double buff1(invSparsity[2]);
	double buff2;

	for (int k = 3; k < J; k++)
	{
		buff2 = buff1;
		buff1 = invSparsity[k];

		invSparsity[k] = (buff1 * (-3.6082 * pow(10, -16)) + buff2 * (0.1716)) - (0.2929 * invSparsity[k-1]) - (0.5858* invSparsity[k-2]) - (0.2929* invSparsity[k-3]);
	}*/
	//peack detection//

	float q(1.0 - (float(gapBetweenFrame) / float(lenthOfFrame)));
	unsigned int h(round(((1.0 - q) * lenthFft)));
	float r(sampleRate / h);
	unsigned int combination_width(ceil(r * lenthFft / sampleRate / 2));

	int alpha(10); // for interval max value between a and b // 10 by default.
	int beta(10); // for interval min value between a and b // 10 by default.
	int a(10);      // for mean max value between a and b // 10 by default.
	int b(10);      // for mean min value between a and b // 10 by default.

	//for a and b :
	//we choose 40 and 40 because each frame is one ms and we don't want to predict outside the 40ms of the middle of the 120ms
	//for alpha and beta :
	//we choose 30 and 30 to let 30ms beetween each detection which allow a tempo of 33 bip par seconde soit 2000 BPM

	float delta(thresholdValue / 100.0);          // threshold for the mean.
	unsigned int p(0);                            // index
	float sum;
	unsigned int beginLoop = std::max(a, alpha);
	unsigned int endLoop = std::max(b, beta);

	//normalization//
	float maxx = (*std::max_element((&invSparsity[0]), (&invSparsity[0] + invSparsity.size())));
	for (int i = 0; i < int(invSparsity.size()); i++)
	{
		invSparsity[i] = invSparsity[i] / maxx;
	}
	//end normalization//

	//peack detection computation//
	for (int i = 40; i < 80; i++) //40ms beacause each frame is 1ms
	{

		sum = 0;
		for (float j = -a;j < b;j++)
		{
			sum += invSparsity[int(double(i) + double(j))];
		}

		if ((invSparsity[i] == float(*std::max_element((invSparsity.begin() + i - alpha), (invSparsity.begin() + i + beta)))) && (invSparsity[i] >= delta + (1.0 / (double(a) + double(b) + 1.0) * sum)) && ((i - p) > combination_width))
		{
			p = i;
			timeOfOnset.push_back(-1);
		}
		else
		{
			timeOfOnset.push_back(invSparsity[i]);
		}
	}
}













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

void Hanning(int lenth, std::vector<float>& hann)
{
	//Hanning Window
	for (int i = 0; i < lenth; i++) {
		double multiplier = 0.5 * (1 - cos(2 * M_PI * i / (lenth - 1.0)));
		hann[i] = multiplier;
	}
	//
}
