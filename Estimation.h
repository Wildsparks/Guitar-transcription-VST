#pragma once
#ifndef ESTIMATION_INCLUDED
#define ESTIMATION_INCLUDED

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Source\package eigen\Eigen\Dense"
#include <random>
#include <iostream>
#include <random>
#include <fstream>
#include <string>

#define NOMBRE_ESTIMATION 500 //500 before opti
struct Model
{
	Eigen::MatrixXd Mu;//2*132
	Eigen::VectorXd W;//1.132
	Eigen::MatrixXd Sigma[132];//2*2*132
};

double JClassifier(Eigen::VectorXd mu, Eigen::MatrixXd C, double P, Eigen::VectorXd X);
Model Calcul_parametre(Eigen::MatrixXd(&featureMatrix)[NOMBRE_ESTIMATION],
	double val1, double val2, double val3, double val4, double val5,
	double val6, double val7, double val8, double val9, double val10,
	double val11, double val12, double val13, double val14, double val15,
	double val16, double val17, double val18, double val19, double val20,
	double val21, double val22, std::default_random_engine generator, std::normal_distribution<double> distribution);

#endif // ESTIMATION_INCLUDED