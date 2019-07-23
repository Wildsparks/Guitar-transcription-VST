#pragma once
#ifndef ONSET_INCLUDED
#define ONSET_INCLUDED

#include "fftw3.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Functions.h"

void Onset(const double* channelData, int bufsize, double sampleRate, std::vector<bool>& timeOfOnset, std::vector<double>& storageActual, std::vector<double>& storagePast, int const& thresholdValue, std::vector<double>& onsetScope, bool& detectionDone, std::vector<double>& storageSparse);

#endif // ONSET_INCLUDED