#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"
#include "fftw3.h"
#define LENGTHFFT 524288//pow(2,19) //19
#define REAL 0
#define IMAG 1//pow(2,19) //19
void Preprocessingold(std::vector<float>& segment)
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
	std::vector<double> Gauss(segment.size());
	Gaussian(segment.size(), Gauss);
	//end - windows//

	for (int i = 0;i < segment.size();++i)
	{
		segmentFreq[i] = segment[i] *(Gauss[i]);
	}

	fft(segmentFreq);

	segmentFreq[0] = 1.0 * segmentFreq[0];
	for (int i = 1; i < segment.size()/2; ++i)
	{
		segmentFreq[i] = 2.0 * segmentFreq[i];
	}
	segmentFreq[segment.size() / 2] = 1.0 * segmentFreq[segment.size() / 2];
	for (int i = segment.size() / 2+1; i < segment.size(); ++i)
	{
		segmentFreq[i] = 0.0;
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
/*
fftw_complex* HilbertFFTForward(std::vector<float>& segment);
double* HilbertFFTBackward(fftw_complex* out);
void LastFFTForward(std::vector<float>& segment, double* hilbertOut);
*/
void Preprocessing(std::vector<float>& segment, std::vector<std::complex<float>>& hilbertOutput)
{
	//=====================================================//
	//==================Hilbert transform==================//
	//=====================================================//

	unsigned int N = segment.size();

	if (N % 2 == 1)
	{
		segment.push_back(0);
		N++;
	}

	fftw_complex* out;
	fftw_complex* in;

	fftw_plan p;
	
	in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N); 
	p = fftw_plan_dft_1d(N, in, out,FFTW_FORWARD, FFTW_ESTIMATE);

	//windows//
	std::vector<double> Gauss(N);
	Gaussian(N, Gauss);
	//end - windows//

	for (int i = 0; i < N; ++i)
	{
		in[i][REAL] = segment[i] * Gauss[i];
		in[i][IMAG] = 0.0;
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	fftw_free(in);

	fftw_complex* hilbertOut;
	hilbertOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	fftw_plan q;
	q = fftw_plan_dft_1d(N, out, hilbertOut, FFTW_BACKWARD, FFTW_ESTIMATE);

	out[0][REAL] = 1.0 * out[0][REAL];
	out[0][IMAG] = 1.0 * out[0][IMAG];
	for (int i = 1; i < N / 2; ++i)
	{
		out[i][REAL] = 2.0 * out[i][REAL];
		out[i][IMAG] = 2.0 * out[i][IMAG];
	}

	out[N / 2][REAL] = 1.0 * out[N / 2][REAL];
	out[N / 2][IMAG] = 1.0 * out[N / 2][IMAG];

	for (int i = N / 2 + 1; i < N; ++i)
	{
		out[i][REAL] = 0.0;
		out[i][IMAG] = 0.0;
	}

	fftw_execute(q);
	fftw_free(out);
	fftw_destroy_plan(q);

	//=====================================================//
	//===============Hilbert transform end=================////NB : FFTW don't normalize fft data : need to be divide by length
	//=====================================================//
	hilbertOutput.clear();
	std::complex<double> valueHilbertOutput(0,0);
	for (int i = 0; i < N; ++i)
	{
		valueHilbertOutput.real(hilbertOut[i][REAL] / double(N));
		valueHilbertOutput.imag(hilbertOut[i][IMAG] / double(N));
		hilbertOutput.push_back(valueHilbertOutput);
	}

	//=====================================================//
	//================Zero padding & FFT===================//
	//=====================================================//

	fftw_complex* in2;
	in2    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * LENGTHFFT);
	fftw_complex* fftOut;
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * LENGTHFFT);
	
	p = fftw_plan_dft_1d(LENGTHFFT, in2, fftOut,FFTW_FORWARD, FFTW_ESTIMATE);

	for (int i = 0; i < N; ++i)
	{
		in2[i][REAL] = hilbertOut[i][REAL] / double(N);
		in2[i][IMAG] = hilbertOut[i][IMAG] / double(N);

	}
	fftw_free(hilbertOut);
	for (int i = N; i < LENGTHFFT; ++i)
	{
		in2[i][REAL] = 0.0;
		in2[i][IMAG] = 0.0;
	}

	fftw_execute(p);
	fftw_free(in2);

	segment.resize(LENGTHFFT / 2);

	for (int i = 0; i < LENGTHFFT / 2; ++i)
	{
		segment[i] = (2.0 / (double)LENGTHFFT * (pow(fftOut[i][IMAG], 2.0) + pow(fftOut[i][REAL], 2.0))); //10.0 * log10() if you want db and *2 can be remove
	}
	fftw_destroy_plan(p);
	fftw_free(fftOut);

	//=====================================================//
	//==============Zero padding & FFT End=================//
	//=====================================================//
	fftw_cleanup();
}
/*
fftw_complex* HilbertFFTForward(std::vector<float>& segment)
{
	fftw_complex in[N], out[N];
	fftw_plan p;
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//windows//
	std::vector<float> Gauss(segment.size());
	Gaussian(segment.size(), Gauss);
	//end - windows//

	for (int i = 0; i < segment.size(); ++i)
	{
		in[i][0] = segment[i] * double(Gauss[i]);
		in[i][1] = 0;
	}
	for (int i = segment.size(); i < N; ++i) //pad
	{
		in[i][0] = 0;
		in[i][1] = 0;
	}

	fftw_execute(p);
	return(out);
	fftw_destroy_plan(p);

}

double* HilbertFFTBackward(fftw_complex* out)
{
	double hilbertOut[N];
	fftw_plan q;
	q = fftw_plan_dft_c2r_1d(N, out, hilbertOut, FFTW_ESTIMATE);

	out[0][0] = 1.0 * out[0][0];
	out[0][1] = 1.0 * out[0][1];
	for (int i = 1; i < N / 2; ++i)
	{
		out[i][0] = 2.0 * out[i][0];
		out[i][1] = 2.0 * out[i][1];
	}

	out[N / 2][0] = 1.0 * out[N / 2][0];
	out[N / 2][1] = 1.0 * out[N / 2][1];

	for (int i = N / 2 + 1; i < N; ++i)
	{
		out[i][0] = 0.0;
		out[i][1] = 0.0;
	}

	fftw_execute(q);
	return(hilbertOut);
	fftw_destroy_plan(q);

}

void LastFFTForward(std::vector<float>& segment, double* hilbertOut)
{
	
	double fftOut[LENGTHFFT];
	fftw_complex in2[LENGTHFFT];
	fftw_plan p;
	p = fftw_plan_dft_c2r_1d(LENGTHFFT, in2, fftOut, FFTW_ESTIMATE);

	for (int i = 0; i < N; ++i)
	{
		in2[i][0] = hilbertOut[i];
		in2[i][1] = 0.0;
	}

	for (int i = N; i < LENGTHFFT; ++i)
	{
		in2[i][0] = 0.0;
		in2[i][1] = 0.0;
	}

	fftw_execute(p);

	segment.resize(LENGTHFFT / 2);

	for (int i = 0; i < LENGTHFFT / 2; ++i)
	{
		segment[i] = (2 * pow(std::abs(fftOut[i]), 2)); //10.0 * log10() if you want db
	}
	fftw_destroy_plan(p);

}
*/