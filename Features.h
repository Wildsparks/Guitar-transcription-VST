#pragma once

#ifndef FEATURES_INCLUDED
#define FEATURES_INCLUDED

#include <iostream>
#include <D:/ecler/Documents/Cours/Ingenieur_4A/Stage/Jacode_III/Source/package eigen/Eigen/Dense>
#include <random>

void featuresExtraction(Eigen::MatrixXd& w0, Eigen::MatrixXd& B, float* channelData, int bufsize);

#endif // FEATURES_INCLUDED
