#pragma once
#ifndef HARMONICSUMMATION_INCLUDED
#define HARMONICSUMMATION_INCLUDED

#include <vector>

double HarmonicSummation(const std::vector<float>& SegmentFFT, unsigned int f0LimitsInf, unsigned int f0LimitsSup, unsigned int nbOfHarmonicsInit, unsigned int sampleRate);

#endif // HARMONICSUMMATION_INCLUDED