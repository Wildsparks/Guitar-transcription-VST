#pragma once
#ifndef ONSET_INCLUDED
#define ONSET_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include "../JuceLibraryCode/JuceHeader.h"

void Onset(const float* channelData, int bufsize, int sampleRate, std::vector<double>& time, int& counter, std::vector<float>& storageActual, std::vector<float>& storagePast);

#endif // ONSET_INCLUDED