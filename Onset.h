#pragma once
#ifndef ONSET_INCLUDED
#define ONSET_INCLUDED

#include "fftw3.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <functional>

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Functions.h"

void Onset(const float* channelData, int bufsize, int sampleRate, std::vector<bool>& timeOfOnset, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast, int const& thresholdValue, std::vector<float>& onsetScope, bool& detectionDone, std::vector<float>& storageSparse);

#endif // ONSET_INCLUDED