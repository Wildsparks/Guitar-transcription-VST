#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"
#include "fftw3.h"
//#define LENGTHFFT (524288/4)//pow(2,19) //19
#define REAL 0
#define IMAG 1

void Preprocessing(std::vector<double>& segment, std::vector<std::complex<double>>& hilbertOutput, int lengthFft)
{
	//=====================================================//
	//==================Hilbert transform==================//
	//=====================================================//

	int N = static_cast<int>(segment.size());

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
	in2    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lengthFft);
	fftw_complex* fftOut;
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lengthFft);
	
	p = fftw_plan_dft_1d(lengthFft, in2, fftOut,FFTW_FORWARD, FFTW_ESTIMATE);

	for (int i = 0; i < N; ++i)
	{
		in2[i][REAL] = hilbertOut[i][REAL] / double(N);
		in2[i][IMAG] = hilbertOut[i][IMAG] / double(N);

	}
	fftw_free(hilbertOut);
	for (int i = N; i < lengthFft; ++i)
	{
		in2[i][REAL] = 0.0;
		in2[i][IMAG] = 0.0;
	}

	fftw_execute(p);
	fftw_free(in2);

	segment.resize(lengthFft / 2);

	for (int i = 0; i < lengthFft / 2; ++i)
	{
		segment[i] = (2.0 / (double)lengthFft * (pow(fftOut[i][IMAG], 2.0) + pow(fftOut[i][REAL], 2.0))); //10.0 * log10() if you want db and *2 can be remove
	}
	fftw_destroy_plan(p);
	fftw_free(fftOut);

	//=====================================================//
	//==============Zero padding & FFT End=================//
	//=====================================================//
	fftw_cleanup();
}
