#pragma once
#ifndef INHARMONIC_INCLUDED
#define INHARMONIC_INCLUDED

#include <vector>

void InharmonicSummation(const std::vector<float>& SegmentFFT, float pitchInitial, unsigned int highNbOfHarmonics, unsigned int sampleRate, double maxBetaGrid, double minBetaGrid, double betaRes, double lengthFFT, float& pitchEstimate, float& BEstimate, float& costFunctionMaxVal);

#endif // INHARMONIC_INCLUDED