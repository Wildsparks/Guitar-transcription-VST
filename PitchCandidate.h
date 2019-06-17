#pragma once
#ifndef PITCHCANDIDATE_INCLUDED
#define PITCHCANDIDATE_INCLUDED

#include <vector>

void PitchCandidate(std::vector<int>& fretCandidate, std::vector<int>& stringCandidate, float observedPitch, std::vector<float> pitchReference, int nbString);

#endif // PITCHCANDIDATE_INCLUDED