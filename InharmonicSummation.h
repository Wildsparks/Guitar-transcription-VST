#pragma once
#ifndef INHARMONIC_INCLUDED
#define INHARMONIC_INCLUDED

#include <vector>

void InharmonicSummation(const std::vector<double>& SegmentFFT, double pitchInitial, int highNbOfHarmonics, double sampleRate, double maxBetaGrid, double minBetaGrid, double betaRes, double lengthFFT, double& pitchEstimate, double& BEstimate, double& costFunctionMaxVal);

#endif // INHARMONIC_INCLUDED