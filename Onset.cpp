#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"


#define TIMEBUFFER 0.04 //40ms of observation

#define REAL 0
#define IMAG 1

void onsetDetector(std::vector<float> const& storage, std::vector<bool> & timeOfOnset, int sampleRate, int const& thresholdValue, std::vector<float> & onsetScope, std::vector<float>& storageSparse);
void onsetDetectorLowCPU(std::vector<float> const& storage, std::vector<bool>& timeOfOnset, int sampleRate, int const& thresholdValue, std::vector<float>& onsetScope);

void Onset(const float* channelData, int bufsize, int sampleRate, std::vector<bool>& timeOfOnset, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast, int const& thresholdValue, std::vector<float>& onsetScope, bool& detectionDone, std::vector<float>& storageSparse)
{

	float nbEchentillon((TIMEBUFFER * sampleRate));        

	//set the first vector of storage//
	if ((storagePast.size() != size_t((double)nbEchentillon * 3.0)) && !detectionDone) //120 milisecond
	{
	
		storagePast.clear();

		while (storagePast.size() != size_t((double)nbEchentillon * 3.0))
		{
			storagePast.push_back(0);
		}

	}
	//done//

	if (detectionDone)
	{
		std::rotate(storagePast.begin(), storagePast.begin() + nbEchentillon, storagePast.end());

		for (int i = 0; i < (nbEchentillon); i++)
		{
			storagePast.pop_back(); //delete last 40ms which are copy
		}

		detectionDone = false;
	}

	for (int i = 0; i < bufsize; i++)
	{

		if (storageActual.size() == nbEchentillon)
		{
	
			for (int j = 0; j < (nbEchentillon); j++)
			{
				storagePast.push_back(storageActual[j]); //add the new 40ms to the old 80ms 
			}

			//120ms analysis from 20ms to 60 ms to get 40ms after and pluck excitation
			onsetDetector(storagePast, timeOfOnset, sampleRate, thresholdValue, onsetScope, storageSparse);

			storageActual.clear();
			storageActual.push_back(channelData[i]);
			detectionDone = true;

		}
		else
		{
			storageActual.push_back(channelData[i]); //40ms
		}
	}

}

void onsetDetector(std::vector<float> const& storage, std::vector<bool>& timeOfOnset, int sampleRate, int const& thresholdValue, std::vector<float>& onsetScope, std::vector<float>& storageSparse)
{

	//time vector contain : 160ms of audio : 4*40ms
	//0.0:0.25 : old 40ms
	//0.25:0.5 : sample to treat
	//0.5:0.75 : futur sample
	//0.75:1.0 : sample for fft length

	int gapBetweenFrame(1);     // in ms
	int lenthOfFrame   (40);    // in ms
	double       gapInSample    (double(gapBetweenFrame) * 1e-3 * double(sampleRate));
	int nbFrame        (ceil(double(storage.size()) / double(gapInSample)));
	int lenthFft       (lenthOfFrame * 1e-3 * sampleRate);
	int beginFrame     (0.5* nbFrame); // future sample to fft
	int stopFrame      (0.25 * nbFrame);

	//windows//
	std::vector<double> hann(lenthFft);
	Hanning(lenthFft, hann);
	//end - windows//

	std::vector <double> frameAbs(lenthFft);
	std::vector <std::vector <double>> spectrogram;

	fftw_complex* frameOut;
	fftw_complex* frameIn;
	frameIn  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lenthFft);
	frameOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lenthFft);

	fftw_plan p;

	p = fftw_plan_dft_1d(lenthFft, frameIn, frameOut, FFTW_FORWARD, FFTW_ESTIMATE);

	//Spectrogram//

	for (int i = beginFrame; i < (nbFrame - stopFrame); ++i)
	{

		for (int j = 0; j < lenthFft; ++j)
		{
			frameIn[j][REAL] = (storage[(i * (double)gapInSample) + j] * double(hann[j]));
			frameIn[j][IMAG] = (0);
		}

		fftw_execute(p);

		for (int k = 0; k < lenthFft; ++k)
		{
			
			frameAbs[k] = std::abs(std::complex<long double>(frameOut[k][REAL], frameOut[k][IMAG]));

		}
		
		spectrogram.push_back(frameAbs);
	}

	fftw_destroy_plan(p);
	fftw_free(frameIn);
	fftw_free(frameOut);
	fftw_cleanup();

	//end - Spectrogram//

	double gamma(0.94); // as in the paper in %
	double J((gamma * double(lenthFft)));
	std::vector <double> invSparsity(0);

	double sumFreq2(0);
	double sumFreq4(0);
	double maxValueOfFrameHarmo(0);

	//inverse sparsity//

	for (int i = 0; i < (spectrogram.size()); ++i)
	{
		maxValueOfFrameHarmo = 0;

		for (int j = 0; j < J; ++j)
		{
			if (maxValueOfFrameHarmo < spectrogram[i][j])
			{
				maxValueOfFrameHarmo = spectrogram[i][j];
			}
		}

		sumFreq2 = 0;
		sumFreq4 = 0;
		for (int j = 0; j < J; ++j)
		{
			spectrogram[i][j] = spectrogram[i][j] / maxValueOfFrameHarmo; //normalize

			if(spectrogram[i][j])

			spectrogram[i][j] = spectrogram[i][j] * spectrogram[i][j];
			sumFreq2 = sumFreq2 + spectrogram[i][j]; //sum(x^2)
			spectrogram[i][j] = spectrogram[i][j] * spectrogram[i][j];
			sumFreq4 = sumFreq4 + spectrogram[i][j]; //sum(x^4)

		}

		if (sumFreq4 == 0)
		{
			invSparsity.push_back(0);
		}
		else
		{
			invSparsity.push_back(sumFreq2 * maxValueOfFrameHarmo / (pow(sumFreq4, 0.25) * pow(J, 0.25)));
		}

	}
	// this is the sparsity method


	//store these value
	for (int i = 0; i < (invSparsity.size()); ++i)
	{
		storageSparse.push_back(invSparsity[i]);
	}

	std::rotate(storageSparse.begin(), storageSparse.begin() + invSparsity.size(), storageSparse.end());

	for (int i = 0; i < (invSparsity.size()); ++i)
	{
		storageSparse.pop_back();
	}
	//

	//take 80ms before

	for (int i = 0; i < (0.5*nbFrame); ++i)
	{
		invSparsity.push_back(storageSparse[floor(storageSparse.size() - 0.75 * nbFrame)+i]);
	}

	//rotate to make it causal

	std::rotate(invSparsity.begin(), invSparsity.begin() + 0.25 * nbFrame, invSparsity.end());


	//peack detection//

	double q(1.0 - (float(gapBetweenFrame) / float(lenthOfFrame)));
	double h(round(((1.0 - q) * lenthFft)));
	double r(sampleRate / h);
	double combination_width(ceil(r * lenthFft / sampleRate / 2.0));

	double alpha(20.0 / (double)gapBetweenFrame); // for interval max value between a and b // 10 by default.
	double beta(20.0 / (double)gapBetweenFrame); // for interval min value between a and b // 10 by default.
	double a(40.0 / (double)gapBetweenFrame);      // for mean max value between a and b // 10 by default.
	double b(40.0 / (double)gapBetweenFrame);      // for mean min value between a and b // 10 by default.

	//for a and b :
	//we choose 40 and 40 because each frame is one ms and we don't want to predict outside the 40ms of the middle of the 120ms
	//for alpha and beta :
	//we choose 30 and 30 to let 30ms beetween each detection which allow a tempo of 33 bip par seconde soit 2000 BPM

	double delta(double(thresholdValue) / 100.0);          // threshold for the mean. 10 advice
	double pointer(0);                            // index
	double sum;
	int beginLoop = std::max(a, alpha);
	int endLoop = std::max(b, beta);

	//normalization//
	double maxx = 0;

	for (int m = 0; m < storageSparse.size(); ++m)
	{
		if (maxx < invSparsity[m])
		{
			maxx = invSparsity[m];
		}
	}

	if (maxx != 0)
	{
		for (int i = 0; i < invSparsity.size(); ++i)
		{
			invSparsity[i] = invSparsity[i] / maxx;
			
		}
	}

	//end normalization//

	timeOfOnset.clear();
	onsetScope.clear();
	
	//peack detection computation//
	for (int i = floor((invSparsity.size()) / 3.0); i < floor((invSparsity.size()) * 2.0 / 3.0); ++i) //40ms beacause each frame is 1ms
	{

		sum = 0;
		for (int j = i - a; j < i + b + 1.0; j++)
		{
			sum = sum + invSparsity[j] + delta;
		}

		maxx = 0;
		for (int m = i - alpha; m <= i + beta; ++m)
		{
			if (maxx < invSparsity[m])
			{
				maxx = invSparsity[m];
			}
		}

		if ((invSparsity[i] == maxx) && (invSparsity[i] >= (sum / (a + b + 1.0 ))) && (((double)i - pointer) > combination_width))
		{
			pointer = i;
			timeOfOnset.push_back(true);
			onsetScope.push_back(-1.0);
		}
		else
		{
			timeOfOnset.push_back(false);
			onsetScope.push_back(invSparsity[i]);
		}

	}


}

void onsetDetectorLowCPU(std::vector<float> const& storage, std::vector<bool>& timeOfOnset, int sampleRate, int const& thresholdValue, std::vector<float>& onsetScope)
{
	//easy method
	unsigned int lengthSegment(storage.size());
	std::vector<float> processSegment(lengthSegment);
	std::vector<float> integrateSegment(lengthSegment);
	std::vector<float> meanSegment(lengthSegment);
	std::vector<float> derivateLogSegment(lengthSegment);

	for (int i = lengthSegment * 0.25 - 150; i < lengthSegment * 0.5 + 151; ++i)
	{
		processSegment[i] = std::abs(storage[i]);
	}

	//infinity norm

	for (int i = lengthSegment * 0.25 - 100; i < lengthSegment * 0.5 + 101; ++i)
	{
		integrateSegment[i] = *std::max_element(processSegment.begin() + i - 50, processSegment.begin() + i + 50);
	}


	//low pass filter

	for (int i = lengthSegment * 0.25 - 50; i < lengthSegment * 0.5 + 51; ++i)
	{
		meanSegment[i] = std::accumulate(integrateSegment.begin() + i - 50, integrateSegment.begin() + i + 50, 0.0);
	}


	//peak creation

	for (int i = lengthSegment * 0.25 - 50; i < lengthSegment * 0.5 + 50; ++i)
	{
		derivateLogSegment[i] = std::abs(meanSegment[i + 1.0] - meanSegment[i]);
	}

	//peack detection//

	unsigned int combination_width(floor(1e-3 * sampleRate));// min 1ms between detection 
	float delta(thresholdValue / 100.0);          // threshold for the mean.
	unsigned int p(0);                            // index

	//normalization//
	float maxx = std::numeric_limits<float>::min();

	for (int i = lengthSegment * 0.25 - 50; i < lengthSegment * 0.5 + 50; ++i)
	{
		if (maxx < derivateLogSegment[i])
		{
			maxx = derivateLogSegment[i];
		}
	}

	for (int i = lengthSegment * 0.25 - 50; i < lengthSegment * 0.5 + 50; ++i)
	{
		derivateLogSegment[i] = (derivateLogSegment[i] / maxx);
	}

	//end normalization//

	timeOfOnset.clear();
	onsetScope.clear();

	//peack detection computation//
	for (int i = floor(0.25 * (double)lengthSegment); i < floor(0.5 * (double)lengthSegment); ++i) //40ms beacause each frame is 1ms
	{

		if ((derivateLogSegment[i] >= delta + std::accumulate(derivateLogSegment.begin() + i - 50, derivateLogSegment.begin() + i + 50, 0.0) / 101) && ((i - p) > combination_width))
		{
			p = i;
			timeOfOnset.push_back(true);
			onsetScope.push_back(derivateLogSegment[i]);
		}
		else
		{
			timeOfOnset.push_back(false);
			onsetScope.push_back(derivateLogSegment[i]);
		}

	}

}

void Onsetold(const float* channelData, int bufsize, int sampleRate, std::vector<bool>& timeOfOnset, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast, int const& thresholdValue, std::vector<float>& onsetScope, bool& detectionDone, std::vector<float>& storageSparse)
{
	float nbEchentillon((TIMEBUFFER * sampleRate));
	int nbBuffers(ceil(nbEchentillon / bufsize));                // number of buffer to store to get TIMEBUFFER miliseconde (ceil get the int sup)

	//set the first vector of storage//
	if ((storagePast.size() != int(double(bufsize) * nbBuffers * 3)) && !detectionDone) //120 milisecond
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
		onsetDetector(storagePast, timeOfOnset, sampleRate, thresholdValue, onsetScope, storageSparse);

		for (int i = 0; i < (bufsize * nbBuffers * 3); i++)
		{
			storagePast[i] = (storagePast[i + double(bufsize) * nbBuffers]); //40ms first milisecond to trash by decay
		}

		for (int i = 0; i < (bufsize * nbBuffers); i++)
		{
			storagePast.pop_back(); //delete last 40ms which are copy
		}

		storageActual.clear();
		for (int i = 0; i < bufsize; i++)
		{
			storageActual.push_back(channelData[i]); //40ms
		}
		detectionDone = true;
		counter++;
	}
	else
	{
		if (detectionDone)
		{
			detectionDone = false;
		}
		for (int i = 0; i < bufsize; i++)
		{
			storageActual.push_back(channelData[i]); //40ms
		}
		counter++;
	}
}