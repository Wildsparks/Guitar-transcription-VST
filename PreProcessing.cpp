#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"
#define LENGTHFFT pow(2,15) //19

void Preprocessing(std::vector<float>& segment)
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
	Hanning(segment.size(), hann);
	//end - windows//

	for (int i = 0;i < segment.size();++i)
	{
		segmentFreq[i] = segment[i]* double(hann[i]);
	}

	fft(segmentFreq);

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

	//=====================================================//
	//================Zero padding & FFT===================//
	//=====================================================//

	while (segment.size() < LENGTHFFT)
	{
		segment.push_back(0);
	}

	segmentFreq.resize(LENGTHFFT);

	for (int i = 0;i < LENGTHFFT;++i)
	{
		segmentFreq[i] = segment[i];
	}

	fft(segmentFreq);

	for (int i = 0;i < LENGTHFFT /2;++i)
	{
		segment[i] =(2 * pow(std::abs(segmentFreq[i]),2)); //10.0 * log10() if you want db
	}
	for (int i = LENGTHFFT / 2;i < LENGTHFFT;++i)
	{
		segment.pop_back();
	}

	//=====================================================//
	//==============Zero padding & FFT End=================//
	//=====================================================//
}
