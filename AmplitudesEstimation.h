#pragma once
#ifndef AMPLITUDESESTIMATION_INCLUDED
#define AMPLITUDESESTIMATION_INCLUDED

#include <D:/ecler/Documents/Cours/Ingenieur_4A/Stage/Jacode_III/Source/package eigen/Eigen/Dense>
#include <vector>

std::vector<double>  AmplitudesEstimation(std::vector<std::complex<float>> const& segment, float f0Features, float betaFeatures, int lengthOfSegment, int sampleRate, int highNbOfHarmonics);
#endif // AMPLITUDESESTIMATION_INCLUDED