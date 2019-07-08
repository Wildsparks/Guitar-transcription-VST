#pragma once
#ifndef ESTIMATION_INCLUDED
#define ESTIMATION_INCLUDED

#include <iostream>
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Source\package eigen\Eigen\Dense"
#include <random>
#define NOMBRE_ESTIMATION 500
struct Model
{
	Eigen::MatrixXd Mu;//2*78
	Eigen::VectorXd W;//1.78
	Eigen::MatrixXd Sigma[78];//2*2*78
};

double JClassifier(Eigen::VectorXd mu, Eigen::MatrixXd C, double P, Eigen::VectorXd X);
Model Calcul_parametre(Eigen::MatrixXd (&featureMatrix)[NOMBRE_ESTIMATION]);

#endif // ESTIMATION_INCLUDED