#pragma once
#ifndef FUNCTIONS_INCLUDED
#define FUNCTIONS_INCLUDED
#include <complex>
#include <valarray>
#include <algorithm> 

#include <vector>

#define M_PI 3.141592653589793238460

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray& x);
void Hanning(int lenth, std::vector<double>& hann);
void ifft(CArray& x);
void Gaussian(int lenth, std::vector<double>& Gaussian);

#endif // FUNCTIONS_INCLUDED