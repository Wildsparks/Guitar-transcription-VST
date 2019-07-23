#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"


#define TIMEBUFFER 0.04 //40ms of observation

#define REAL 0
#define IMAG 1

void onsetDetector(std::vector<double> const& storage, std::vector<bool> & timeOfOnset, double sampleRate, int const& thresholdValue, std::vector<double> & onsetScope, std::vector<double> & storageSparse);
void onsetDetectorLowCPU(std::vector<double> const& storage, std::vector<bool>& timeOfOnset, double sampleRate, int const& thresholdValue, std::vector<double>& onsetScope, std::vector<double>& storageSparse);

void Onset(const double* channelData, int bufsize, double sampleRate, std::vector<bool>& timeOfOnset, std::vector<double>& storageActual, std::vector<double>& storagePast, int const& thresholdValue, std::vector<double>& onsetScope, bool& detectionDone, std::vector<double>& storageSparse)
{

	double nbEchentillon((TIMEBUFFER * sampleRate));

	//set the first vector of storage//
	if ((storagePast.size() != size_t(nbEchentillon * 3.0)) && !detectionDone) //120 milisecond
	{

		storagePast.clear();

		while (storagePast.size() != size_t(nbEchentillon * 3.0))
		{
			storagePast.push_back(0);
		}

	}
	//done//

	if (detectionDone)
	{
		std::rotate(storagePast.begin(),storagePast.begin() + static_cast<__int64>(nbEchentillon), storagePast.end());

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

void onsetDetector(std::vector<double> const& storage, std::vector<bool>& timeOfOnset, double sampleRate, int const& thresholdValue, std::vector<double>& onsetScope, std::vector<double>& storageSparse)
{

	//time vector contain : 160ms of audio : 4*40ms
	//0.0:0.25 : old 40ms
	//0.25:0.5 : sample to treat
	//0.5:0.75 : futur sample
	//0.75:1.0 : sample for fft length

	int     gapBetweenFrame(1);     // in ms
	int     lenthOfFrame   (40);    // in ms
	double  gapInSample    (gapBetweenFrame * 1e-3 * sampleRate);
	double  nbFrame        ((storage.size() / gapInSample));
	int     lenthFft       (static_cast<int>(lenthOfFrame * 1e-3 * sampleRate));
	int     beginFrame     (static_cast<int>(0.5 * nbFrame)); // future sample to fft
	int     stopFrame      (static_cast<int>(0.25 * nbFrame));

	//windows//
	std::vector<double> hann(lenthFft);
	Hanning(lenthFft, hann);
	//end - windows//

	std::vector <double> frameAbs(lenthFft);
	std::vector <std::vector <double>> spectrogram;

	fftw_complex* frameOut;
	fftw_complex* frameIn;
	frameIn = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lenthFft * 40);
	frameOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lenthFft * 40);

	fftw_plan p;

	//p = fftw_plan_dft_1d(lenthFft, frameIn, frameOut, FFTW_FORWARD, FFTW_ESTIMATE);


	int rank = 1;
	int n[] = { lenthFft };
	int howmany = 40;
	int idist = lenthFft;
	int odist = lenthFft;      /*  the distance in memory between the first element of the first array and the first element of the second array */
	int istride = 1;           /* array is contiguous in memory */
	int ostride = 1;
	int* inembed = n, * onembed = n;

	p = fftw_plan_many_dft( rank, n, howmany,
		                    frameIn, inembed,
                            istride, idist,
                            frameOut,onembed,
                            ostride, odist,
							FFTW_FORWARD, FFTW_ESTIMATE);



	//Spectrogram//

	for (int i = beginFrame; i < (nbFrame - stopFrame); ++i)
	{
		for (int j = 0; j < lenthFft; ++j)
		{
			frameIn[j+ lenthFft * (i - beginFrame)][REAL] = (storage[static_cast<__int64>((i * gapInSample) + j)] * (hann[j]));
			frameIn[j+ lenthFft * (i - beginFrame)][IMAG] = (0);
		}
	}

	fftw_execute(p);

	for (int i = 0; i < (nbFrame - stopFrame - beginFrame); ++i)
	{
		for (int k = 0; k < lenthFft; ++k)
		{

			frameAbs[k] = std::abs(std::complex<long double>(frameOut[k + lenthFft * i][REAL], frameOut[k + lenthFft * i][IMAG]));

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

	double sumFreq2(0.0);
	double sumFreq4(0.0);
	double maxValueOfFrameHarmo(0.0);

	//inverse sparsity//

	for (int i = 0; i < (spectrogram.size()); ++i)
	{
		maxValueOfFrameHarmo = 0.0;

		for (int j = 0; j < J; ++j)
		{
			if (maxValueOfFrameHarmo < spectrogram[i][j])
			{
				maxValueOfFrameHarmo = spectrogram[i][j];
			}
		}

		sumFreq2 = 0.0;
		sumFreq4 = 0.0;
		for (int j = 0; j < J; ++j)
		{
			spectrogram[i][j] = spectrogram[i][j] / maxValueOfFrameHarmo; //normalize

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

	for (int i = 0; i < (0.5 * nbFrame); ++i)
	{
		invSparsity.push_back(storageSparse[storageSparse.size() - static_cast<int>(0.75 * nbFrame) + i]);
	}

	//rotate to make it causal

	std::rotate(invSparsity.begin(), invSparsity.begin() + static_cast<int>(0.25 * nbFrame), invSparsity.end());


	//peack detection//

	double q(1.0 - (double(gapBetweenFrame) / double(lenthOfFrame)));
	double h(round((1.0 - q) * double(lenthFft)));
	double r(sampleRate / h);
	double combination_width(ceil(r * lenthFft / sampleRate / 2.0));

	double alpha ((20.0 / (double)gapBetweenFrame));      // for interval max value between a and b // 20 by default.
	double beta  ((20.0 / (double)gapBetweenFrame));      // for interval min value between a and b // 20 by default.
	double a     ((40.0 / (double)gapBetweenFrame));      // for mean max value between a and b // 40 by default.
	double b     ((00.0 / (double)gapBetweenFrame));      // for mean min value between a and b // 40 by default.

	double delta(double(thresholdValue) / 100.0);          // threshold for the mean. 10 advice
	int    pointer(0);                                     // index
	double sum;

	//normalization//
	double maxx = 0;

	for (int m = 0; m < storageSparse.size(); ++m)
	{
		if (maxx < storageSparse[m])
		{
			maxx = storageSparse[m];
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
	for (int i = static_cast<int>((invSparsity.size()) / 3.0); i < static_cast<int>((invSparsity.size()) * 2.0 / 3.0); ++i) //40ms beacause each frame is 1ms
	{

		sum = 0;
		for (int j = static_cast<int>(i - a); j < static_cast<int>(i + b + 1.0); j++)
		{
			sum = sum + invSparsity[j] + delta;
		}

		maxx = 0;
		for (int m = static_cast<int>(i - alpha); m <= static_cast<int>(i + beta); ++m)
		{
			if (maxx < invSparsity[m])
			{
				maxx = invSparsity[m];
			}
		}

		if ((invSparsity[i] == maxx) && (invSparsity[i] >= (sum / (a + b + 1.0))) && (double(i - pointer) > combination_width))
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


void onsetDetectorLowCPU(std::vector<double> const& storage, std::vector<bool>& timeOfOnset, double sampleRate, int const& thresholdValue, std::vector<double>& onsetScope, std::vector<double>& storageSparse)
{

	//easy method
	int gap = 10; //ms

	int lengthSegment = static_cast<int>( storage.size() );

	double scale = (gap * 0.001 * sampleRate);
	int nbFrame = static_cast<int>(lengthSegment / scale);
	int lengthOfAnalyse = static_cast<int>(round(scale));

	std::vector<double> processSegment(lengthSegment);
	std::vector<double> integrateSegment(lengthSegment);
	std::vector<double> meanSegment(lengthSegment);
	std::vector<double> derivateLogSegment(lengthSegment);

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse*3.0); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse*3.0+1.0); ++i)
	{
		processSegment[i] = std::abs(storage[i]);
	}

	//infinity norm

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse*2.0); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse*2.0 +1.0); ++i)
	{
		integrateSegment[i] = *std::max_element(processSegment.begin() + i - lengthOfAnalyse, processSegment.begin() + i + lengthOfAnalyse);
	}


	//low pass filter

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse + 1.0); ++i)
	{
		meanSegment[i] = std::accumulate(integrateSegment.begin() + i - lengthOfAnalyse, integrateSegment.begin() + i + lengthOfAnalyse, 0.0);
	}


	//peak creation

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse); ++i)
	{
		derivateLogSegment[i] = (meanSegment[i + 1] - meanSegment[i]);
		if(derivateLogSegment[i]<0)
		{
			derivateLogSegment[i] = 0;
		}
	}

	//peack detection//

	//normalization//
	
	double sum(0);
	double maxx = std::numeric_limits<float>::min();

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse); ++i)
	{
		if (maxx < derivateLogSegment[i])
		{
			maxx = derivateLogSegment[i];
		}
	}

	//store this value

	storageSparse.push_back(maxx);

	std::rotate(storageSparse.begin(), storageSparse.begin() + 1, storageSparse.end());

	storageSparse.pop_back();

	maxx = std::numeric_limits<float>::min();
	//
	for (int i = 0; i < storageSparse.size(); ++i)
	{
		if (maxx < storageSparse[i])
		{
			maxx = storageSparse[i];
		}
	}
	//

	for (int i = static_cast<int>(lengthSegment * 0.25 - lengthOfAnalyse); i < static_cast<int>(lengthSegment * 0.5 + lengthOfAnalyse); ++i)
	{
		derivateLogSegment[i] = (derivateLogSegment[i] / maxx);
	}

	//end normalization//

	int combination_width(lengthOfAnalyse);// min 1ms between detection 
	double delta(double(0.10));//thresholdValue) / 100.0);           // threshold for the mean.
	double pointer(0);                            // index

	double alpha(lengthOfAnalyse);      // for interval max value between a and b // 20 by default.
	double beta(lengthOfAnalyse);      // for interval min value between a and b // 20 by default.
	double a(lengthOfAnalyse);      // for mean max value between a and b // 40 by default.
	double b(lengthOfAnalyse);      // for mean min value between a and b // 40 by default.

	timeOfOnset.clear();
	onsetScope.clear();

	//peack detection computation//
	for (int i = static_cast<int>(0.25 * (double)lengthSegment); i < static_cast<int>(0.5 * (double)lengthSegment); ++i)
	{

		sum = 0;
		for (int j = static_cast<int>(i - a); j < static_cast<int>(i + b + 1.0); ++j)
		{
			sum = sum + derivateLogSegment[j] + delta;
		}

		maxx = 0;
		for (int m = static_cast<int>(i - alpha); m <= static_cast<int>(i + beta); ++m)
		{
			if (maxx < derivateLogSegment[m])
			{
				maxx = derivateLogSegment[m];
			}
		}

		if ((derivateLogSegment[i] == maxx) && (derivateLogSegment[i] >= (sum / (a + b + 1.0))) && (double(double(i) - pointer) > combination_width))
		{
			pointer = i;
			timeOfOnset.push_back(true);
			onsetScope.push_back(-1.0);
		}
		else
		{
			timeOfOnset.push_back(false);
			onsetScope.push_back(derivateLogSegment[i]);
		}

	}

}