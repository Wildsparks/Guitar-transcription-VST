#pragma once
#ifndef FEATURES_INCLUDED
#define FEATURES_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <random>

void featuresExtraction(Eigen::MatrixXd& w0, Eigen::MatrixXd& B, float* channelData, int bufsize);

#endif // FEATURES_INCLUDED