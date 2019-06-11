#pragma once
#ifndef ONSET_INCLUDED
#define ONSET_INCLUDED

#include <iostream>

#include <complex>
#include <valarray>
#include <algorithm> 

#include <vector>
#include <numeric>
#include <string>
#include <functional>

void Onset(const float* channelData, int bufsize, int sampleRate, std::vector<double>& timeOfOnset, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast, int const& thresholdValue);

#endif // ONSET_INCLUDED
