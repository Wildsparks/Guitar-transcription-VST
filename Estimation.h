#pragma once
#ifndef ESTIMATION_INCLUDED
#define ESTIMATION_INCLUDED

#include <iostream>
#include <Eigen/Dense>
#include <random>

struct Model
{
	Eigen::MatrixXd Mu;//2*78
	Eigen::VectorXd W;//1.78
	Eigen::MatrixXd Sigma[78];//2*2*78
};

double JClassifier(Eigen::VectorXd mu, Eigen::MatrixXd C, double P, Eigen::VectorXd X);
Model Calcul_parametre();

#endif // ESTIMATION_INCLUDED