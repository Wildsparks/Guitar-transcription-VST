#pragma once
#ifndef HARMONICSUMMATION_INCLUDED
#define HARMONICSUMMATION_INCLUDED

#include <vector>

double HarmonicSummation(const std::vector<float>& SegmentFFT, int f0LimitsInf, int f0LimitsSup, int nbOfHarmonicsInit, double sampleRate);

#endif // HARMONICSUMMATION_INCLUDED