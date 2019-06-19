#pragma once
#ifndef INHARMONIC_INCLUDED
#define INHARMONIC_INCLUDED

#include <vector>

void InharmonicSummation(std::vector<float>& const SegmentFFT, float pitchInitial, unsigned int highNbOfHarmonics, unsigned int sampleRate, double maxBetaGrid, double minBetaGrid, float betaRes, float lengthFFT, float& pitchEstimate, float& BEstimate, float& costFunctionMaxVal);

#endif // INHARMONIC_INCLUDED