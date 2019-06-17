#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"
#include <iostream>
#include <random>

using namespace Eigen;

#define NOMBRE_ESTIMATION 500 //5000
# define M_PI 3.14159265358979323846

struct Model;

double std_dev(ArrayXd vec)
{
	double std_dev = std::sqrt((vec - vec.mean()).square().sum() / (vec.size() - 1.0));
	return std_dev;
}

MatrixXd covariance(MatrixXd mat)
{

	MatrixXd centered = mat.rowwise() - mat.colwise().mean();
	MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows() - 1.0);
	return cov;
}

double Vect_Mat_Vect(VectorXd A, MatrixXd C, VectorXd B)
{
	return B(1) * C(1, 1) * A(1) + B(0) * C(1, 0) * A(1) + B(1) * C(0.1) * A(0) + B(0) * C(0, 0) * A(0);
}

double JClassifier(VectorXd mu, MatrixXd C, double P, VectorXd X)
{
	MatrixXd CI = C.inverse();

	return -log(C.determinant()) + 2 * log(P) - Vect_Mat_Vect(mu, CI, mu) + 2 * Vect_Mat_Vect(X, CI, mu) - Vect_Mat_Vect(X, CI, X);
}


Model Calcul_parametre()
{
	//randn
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);
	//-------------------------------//
	double percentageStd = 0.00005;
	double number = distribution(generator);
	double P(1.0 / 3.0);
	double force(0.05);
	double Esteel = 2.27e11; // Young's modulus for steel
	double G = 79.3e9; // known shear modulus steel
	double rhoCore = 7950; //[kg/m^3]
	double rhoWrapping = 6000; //[kg/m^3]

	VectorXd distanceFactor(6);
	distanceFactor << 1, 2, 9, 4, 5, 2;

	VectorXd L0(6);
	for (int i = 0; i < 6; i++)
	{
		L0(i) = 0.6411 + distanceFactor(i) * 1e-3; //distance from nut to bridge
	}

	int numStrings(6);
	VectorXd strNdx(numStrings);
	for (int i = 0; i < numStrings; i++)
	{
		strNdx(i) = i++;
	}

	int numFrets(12);
	VectorXd fretNdx(numFrets + 1);
	for (int i = 0; i <= numFrets; i++)
	{
		fretNdx(i) = i;
	}

	//dFull
	VectorXd thicknessFactor(6);
	thicknessFactor << 0.010, 0.013, 0.0172, 0.026, 0.036, 0.046;

	VectorXd dFull(6);
	dFull = thicknessFactor * 0.0254; //full diameter for the  strings


	VectorXd thicknessFactorCore(6);
	thicknessFactorCore << 0.010, 0.013, 0.0172, 0.0146, 0.016, 0.018;


	VectorXd dCore(6);
	dCore = thicknessFactorCore * 0.0254; //core diameter without wrapping


	VectorXd T0(6);
	T0 << 16.5, 16, 17.5, 18.5, 20, 17;
	T0 *= 4.45;

	////-------------------------------------------////

	MatrixXd BIntrinsic[NOMBRE_ESTIMATION];
	MatrixXd BPluck[NOMBRE_ESTIMATION];
	MatrixXd f0[NOMBRE_ESTIMATION];
	MatrixXd BMixed[NOMBRE_ESTIMATION];

	////-------------------------------------------////

	for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; nombre_estimation++)
	{
		//add noise to length
		MatrixXd mL0(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				mL0(j, i) = L0(j) * pow(2, -fretNdx(i) / 12);
			}
		}

		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				mL0(j, i) = mL0(j, i) + percentageStd * distribution(generator) * mL0(j, i);
			}
		}

		//noise to T0
		for (int j = 0; j < numStrings; j++)
		{
			T0(j) = T0(j) + percentageStd * distribution(generator) * T0(j);
		}

		//noise to dCore
		for (int j = 0; j < numStrings; j++)
		{
			dCore(j) = dCore(j) + percentageStd * distribution(generator) * dCore(j);
		}

		VectorXd dWrapping((dFull - dCore) * 0.5);

		//noise to dWrapping
		for (int j = 0; j < numStrings; j++)
		{
			dWrapping(j) = dWrapping(j) + percentageStd * distribution(generator) * dWrapping(j);
		}

		//Add noise to plucking force
		force = force + percentageStd * distribution(generator) * force;


		///--------------------------------------------///
		///--------------CALCUL B & F0-----------------///
		///--------------------------------------------///

		VectorXd ACore(6);
		for (int j = 0; j < numStrings; j++)
		{
			ACore(j) = M_PI * pow((dCore(j) / 2), 2); /*(pi* dCore. ^ 2) / 4;% Crosssection core[m ^ 2]*/
		}

		// Mass - per - unit length


		VectorXd mu(6);
		for (int j = 0; j < numStrings; j++)
		{
			mu(j) = ACore(j) * rhoCore + rhoWrapping * (pow((2 * dWrapping(j) + dCore(j)), 2) - pow(dCore(j), 2)) * (M_PI / 4);
		}

		//Calculating deltaL from material properties

		VectorXd D(dCore + dWrapping);

		VectorXd TcOverTw(6);
		for (int j = 0; j < numStrings; j++)
		{
			TcOverTw(j) = (8.0 * ACore(j) * pow(D(j), 3.0) * Esteel) / (G * pow(dWrapping(j), 5.0));
		}

		VectorXd Tc(6);
		for (int j = 0; j < numStrings; j++)
		{
			Tc(j) = T0(j) / ((1 / TcOverTw(j)) + 1.0);
		}

		VectorXd Tw(6);
		for (int j = 0; j < numStrings; j++)
		{
			Tw(j) = T0(j) / (TcOverTw(j) + 1.0);
		}

		MatrixXd deltaL(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaL(j, i) = (mL0(j, i) * Tc(j)) / (ACore(j) * Esteel + Tc(j));
			}
		}

		//calculate effective Young's modulus for all strings

		MatrixXd E(6, 13);
		MatrixXd Eeff(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				E(j, i) = (Tc(j) / ACore(j)) / (deltaL(j, i) / (mL0(j, i) - deltaL(j, i)));
				Eeff(j, i) = (T0(j) / ACore(j)) / (deltaL(j, i) / (mL0(j, i) - deltaL(j, i)));
			}
		}

		//Transverse Displacement [Abbot, Strings in the 16th and 17th Centuries, 1974, Appendix 3]
		MatrixXd deltaP(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaP(j, i) = pow(((mL0(j, i) * P * force * (1 - P)) / T0(j)), 2);
			}
		}

		//Length extension

		MatrixXd deltaDeltaL(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaDeltaL(j, i) = sqrt(pow((P * mL0(j, i)), 2) + pow(deltaP(j, i), 2)) + sqrt(pow((1 - P) * mL0(j, i), 2) + pow(deltaP(j, i), 2)) - mL0(j, i);
			}
		}

		//Displacement as if plucked in the middle (retaining the length extension)

		MatrixXd deltaHalf(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaHalf(j, i) = sqrt((pow(deltaDeltaL(j, i), 2.0) + deltaDeltaL(j, i) * mL0(j, i) * 2.0) / 4.0);
			}
		}

		//Calculate inharmonicity factor

		MatrixXd K(6, 13);
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				K(j, i) = (pow(M_PI, 3) * Eeff(j, i) * pow(dCore(j), 2)) / (16 * T0(j) * pow(mL0(j, i), 2));
			}
		}
		//K = -(pi^3 * Eeff .* dCore.^2) ./ (16 * T0 .* L0.^1);

		//--------------------------------------------------//

		//B is a sum of intrinsic and pluck deflection
		//	MatrixXd BIntrinsic[NOMBRE_ESTIMATION]  //
		//  MatrixXd BPluck[NOMBRE_ESTIMATION];     //
		//  MatrixXd f0[NOMBRE_ESTIMATION];         //
		//  MatrixXd BMixed[NOMBRE_ESTIMATION];     //

		MatrixXd BIntrinsicForOneEstimation(6, 13);
		MatrixXd BPluckForOneEstimation(6, 13);
		MatrixXd f0ForOneEstimation(6, 13);

		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				BIntrinsicForOneEstimation(j, i) = (K(j, i) / 4) * pow(dCore(j), 2);
				BPluckForOneEstimation(j, i) = ((K(j, i) * 3) / 8) * pow(deltaHalf(j, i), 2);
				f0ForOneEstimation(j, i) = sqrt(T0(j) / mu(j)) / mL0(j, i) / 2;
			}
		}

		BIntrinsic[nombre_estimation] = BIntrinsicForOneEstimation;
		BPluck[nombre_estimation] = BPluckForOneEstimation;
		f0[nombre_estimation] = f0ForOneEstimation;

		BMixed[nombre_estimation] = BIntrinsicForOneEstimation + BPluckForOneEstimation;

	}

	//-----END-LOOP-----//

	MatrixXd BSimulations[NOMBRE_ESTIMATION];
	MatrixXd pitchSimulations[NOMBRE_ESTIMATION];
	MatrixXd w0[NOMBRE_ESTIMATION];

	for (int i = 0; i < NOMBRE_ESTIMATION; i++)
	{
		BSimulations[i] = BMixed[i];
		pitchSimulations[i] = f0[i];
		w0[i] = f0[i];
	}

	//------CLASSIFIER-----//

	MatrixXd featureMatrix[NOMBRE_ESTIMATION];
	Model classifier;

	int kk = 0;

	for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; nombre_estimation++)
	{
		MatrixXd featureMatrixForOneEstimation(2, 78);
		kk = 0;
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = numStrings - 1; j >= 0; j--)
			{

				featureMatrixForOneEstimation(0, kk) = w0[nombre_estimation](j, i);
				featureMatrixForOneEstimation(1, kk) = BMixed[nombre_estimation](j, i);
				kk = kk + 1;

			}
		}
		featureMatrix[nombre_estimation] = featureMatrixForOneEstimation;
	}

	kk = 0;
	MatrixXd meanFeatureMatrix(2, 78);
	for (int i = 0; i <= numFrets; ++i)
	{
		for (int j = numStrings - 1; j >= 0; --j)
		{
			meanFeatureMatrix(0, kk) = 0;
			meanFeatureMatrix(1, kk) = 0;
			for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; ++nombre_estimation)
			{
				meanFeatureMatrix(0, kk) += featureMatrix[nombre_estimation](0, kk);
				meanFeatureMatrix(1, kk) += featureMatrix[nombre_estimation](1, kk);

			}

			meanFeatureMatrix(0, kk) /= NOMBRE_ESTIMATION;
			meanFeatureMatrix(1, kk) /= NOMBRE_ESTIMATION;

			kk = kk + 1;
		}
	}
	classifier.Mu = meanFeatureMatrix;

	kk = 0;
	VectorXd W(78);
	for (int i = 0; i <= numFrets; i++)
	{
		for (int j = numStrings - 1; j >= 0; j--)
		{
			MatrixXd covMatrixTemp(NOMBRE_ESTIMATION, 2);

			for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; nombre_estimation++)
			{
				covMatrixTemp(nombre_estimation, 0) = featureMatrix[nombre_estimation](0, kk);
				covMatrixTemp(nombre_estimation, 1) = featureMatrix[nombre_estimation](1, kk);

			}

			classifier.Sigma[kk] = covariance(covMatrixTemp);
			W(kk) = double(1.0 / (numStrings * (numFrets + 1.0)));

			kk = kk + 1;
		}
	}

	classifier.W = W;

	return classifier;

}
