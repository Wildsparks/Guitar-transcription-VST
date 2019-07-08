#pragma once
#ifndef PITCHCANDIDATE_INCLUDED
#define PITCHCANDIDATE_INCLUDED

#include <vector>

void PitchCandidate(std::vector<int>& finalCandidate, float observedPitch, std::vector<double> pitchReference, int nbString);

#endif // PITCHCANDIDATE_INCLUDED