#pragma once
#ifndef AMPLITUDESESTIMATION_INCLUDED
#define AMPLITUDESESTIMATION_INCLUDED

#include <D:/ecler/Documents/Cours/Ingenieur_4A/Stage/Jacode_III/Source/package eigen/Eigen/Dense>
#include <vector>

std::vector<float>  AmplitudesEstimation(std::vector<float> const& segment, double f0Features, double betaFeatures, int lengthOfSegment, unsigned int sampleRate, unsigned int highNbOfHarmonics);
#endif // AMPLITUDESESTIMATION_INCLUDED