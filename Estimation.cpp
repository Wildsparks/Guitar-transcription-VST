#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"

using namespace Eigen;

# define M_PI 3.1415926535897932

struct Model;


MatrixXd covariance(MatrixXd mat)
{
	MatrixXd centered = mat.rowwise() - mat.colwise().mean();
	MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows() - 1.0);
	return cov;
}

double Vect_Mat_Vect(VectorXd A, MatrixXd C, VectorXd B)
{
	return(A.transpose() * C * B);
}

double JClassifier(VectorXd mu, MatrixXd C, double P, VectorXd X)
{
	MatrixXd CI = C.inverse();

	return -log(C.determinant()) + 2.0 * log(P) - Vect_Mat_Vect(mu, CI, mu) + 2.0 * Vect_Mat_Vect(X, CI, mu) - Vect_Mat_Vect(X, CI, X);
}


Model Calcul_parametre(MatrixXd (&featureMatrix)[NOMBRE_ESTIMATION])
{

	std::string line;
	std::ifstream myfile;
	myfile.open("random_number.txt");

	if (myfile.is_open())
	{}
	else std::cout << "Unable to open file call random_number.text, please join it to the .cpp file folder";

	//randn
	//std::default_random_engine generator; //not same as matlab
	//std::normal_distribution<double> distribution(0.0, 1.0);
	//-------------------------------//

	VectorXd distanceFactor(6);
	VectorXd L0(6);
	int numStrings(6);
	VectorXd strNdx(numStrings);
	int numFrets(12);
	VectorXd fretNdx(numFrets + 1);

	//dFull
	VectorXd thicknessFactor(6);
	VectorXd dFull(6);
	VectorXd thicknessFactorCore(6);
	VectorXd dCore(6);
	VectorXd T0(6);

	////-------------------------------------------////

	MatrixXd f0[NOMBRE_ESTIMATION];
	MatrixXd BMixed[NOMBRE_ESTIMATION];

	////-------------------------------------------////

	MatrixXd mL0(6, 13);
	VectorXd ACore(6);
	VectorXd dWrapping(6);
	VectorXd mu(6);
	VectorXd D(6);
	VectorXd TcOverTw(6);
	VectorXd Tw(6);
	VectorXd Tc(6);
	MatrixXd deltaL(6, 13);
	MatrixXd E(6, 13);
	MatrixXd Eeff(6, 13);
	MatrixXd deltaP(6, 13);
	MatrixXd deltaDeltaL(6, 13);
	MatrixXd deltaHalf(6, 13);
	MatrixXd K(6, 13);
	MatrixXd BIntrinsicForOneEstimation(6, 13);
	MatrixXd BPluckForOneEstimation(6, 13);
	MatrixXd f0ForOneEstimation(6, 13);

	for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; ++nombre_estimation)
	{
		//initialisation
		double percentageStd = 0.005;
		double P(1.0 / 3.0);
		double force(0.05);
		double Esteel = 2.27e11; // Young's modulus for steel
		double G = 79.3e9; // known shear modulus steel
		double rhoCore = 7950.0; //[kg/m^3]
		double rhoWrapping = 6000.0; //[kg/m^3]

		for (int i = 0; i < numStrings; i++)
		{
			strNdx(i) = i++;
		}

		for (int i = 0; i <= numFrets; i++)
		{
			fretNdx(i) = i;
		}

		distanceFactor << 1.0, 2.0, 9.0, 4.0, 5.0, 2.0;

		for (int i = 0; i < 6; i++)
		{
			L0(i) = 0.6411 + distanceFactor(i) * 1e-3; //distance from nut to bridge
		}

		//thicknessFactor << 0.00989, 0.0132, 0.0174, 0.026, 0.036, 0.0459;
		thicknessFactor << 0.010, 0.013, 0.0172, 0.026, 0.036, 0.046;
		dFull = thicknessFactor * 0.0254; //full diameter for the  strings
		thicknessFactorCore << 0.010, 0.013, 0.0172, 0.0146, 0.016, 0.018;
		dCore = thicknessFactorCore * 0.0254; //core diameter without wrapping
		T0 << 16.5, 16.0, 17.5, 18.5, 20, 17;
		T0 *= 4.45;

		//add noise to length
		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				mL0(j, i) = L0(j) * pow(2.0, -(fretNdx(i)) / 12.0);
			}
		}

		for (int j = 0; j < numStrings; j++)
		{
			for (int i = 0; i <= numFrets; i++)
			{
				getline(myfile, line); // for random
				mL0(j, i) = mL0(j, i) + percentageStd * std::stod(line) * mL0(j, i);
			}
		}

		//noise to T0
		for (int j = 0; j < numStrings; j++)
		{
			getline(myfile, line); // for random
			T0(j) = T0(j) + percentageStd * std::stod(line) * T0(j);
		}

		//noise to dCore
		for (int j = 0; j < numStrings; j++)
		{
			getline(myfile, line); // for random
			dCore(j) = dCore(j) + percentageStd * std::stod(line) * dCore(j);
		}

		for (int j = 0; j < numStrings; j++)
		{
			dWrapping(j) = ((dFull(j) - dCore(j)) / 2.0);
		}

		//noise to dWrapping
		for (int j = 0; j < numStrings; j++)
		{
			getline(myfile, line); // for random
			dWrapping(j) = dWrapping(j) + percentageStd * std::stod(line) * dWrapping(j);
		}

		//Add noise to plucking force
		getline(myfile, line); // for random
		force = force + percentageStd * std::stod(line) * force;


		///--------------------------------------------///
		///--------------CALCUL B & F0-----------------///
		///--------------------------------------------///

		
		for (int j = 0; j < numStrings; j++)
		{
			ACore(j) = M_PI * pow((dCore(j) / 2.0), 2.0); /*(pi* dCore. ^ 2) / 4;% Crosssection core[m ^ 2]*/
		}

		// Mass - per - unit length


		for (int j = 0; j < numStrings; j++)
		{
			mu(j) = ACore(j) * rhoCore + rhoWrapping * (pow((2.0 * dWrapping(j) + dCore(j)), 2.0) - pow(dCore(j), 2.0)) * (M_PI / 4.0);
		}

		//Calculating deltaL from material properties

		D=(dCore + dWrapping);

		
		for (int j = 0; j < numStrings; j++)
		{
			TcOverTw(j) = (8.0 * ACore(j) * pow(D(j), 3.0) * Esteel) / (G * pow(dWrapping(j), 5.0));
		}

		
		for (int j = 0; j < numStrings; j++)
		{
			Tc(j) = T0(j) / ((1.0 / TcOverTw(j)) + 1.0);
		}

		
		for (int j = 0; j < numStrings; j++)
		{
			Tw(j) = T0(j) / (TcOverTw(j) + 1.0);
		}

		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaL(j, i) = (mL0(j, i) * Tc(j)) / (ACore(j) * Esteel + Tc(j));
			}
		}

		//calculate effective Young's modulus for all strings

		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				E(j, i) = (Tc(j) / ACore(j)) / (deltaL(j, i) / (mL0(j, i) - deltaL(j, i)));
				Eeff(j, i) = (T0(j) / ACore(j)) / (deltaL(j, i) / (mL0(j, i) - deltaL(j, i)));
			}
		}

		//Transverse Displacement [Abbot, Strings in the 16th and 17th Centuries, 1974, Appendix 3]
		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaP(j, i) = pow(((mL0(j, i) * P * force * (1.0 - P)) / T0(j)), 2.0);
			}
		}

		//Length extension

		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaDeltaL(j, i) = std::sqrt(pow((P * mL0(j, i)), 2.0) + pow(deltaP(j, i), 2.0)) + std::sqrt(pow((1.0 - P) * mL0(j, i), 2.0) + pow(deltaP(j, i), 2.0)) - mL0(j, i);
			}
		}

		//Displacement as if plucked in the middle (retaining the length extension)

		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				deltaHalf(j, i) = sqrt((pow(deltaDeltaL(j, i), 2.0) + deltaDeltaL(j, i) * mL0(j, i) * 2.0) / 4.0);
			}
		}

		//Calculate inharmonicity factor

		
		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				K(j, i) = (pow(M_PI, 3.0) * Eeff(j, i) * pow(dCore(j), 2.0)) / (16.0 * T0(j) * pow(mL0(j, i), 2.0));
			}
		}
		//K = -(pi^3 * Eeff .* dCore.^2) ./ (16 * T0 .* L0.^1);

		//--------------------------------------------------//

		//B is a sum of intrinsic and pluck deflection
		//	MatrixXd BIntrinsic[NOMBRE_ESTIMATION]  //
		//  MatrixXd BPluck[NOMBRE_ESTIMATION];     //
		//  MatrixXd f0[NOMBRE_ESTIMATION];         //
		//  MatrixXd BMixed[NOMBRE_ESTIMATION];     //

		for (int i = 0; i <= numFrets; i++)
		{
			for (int j = 0; j < numStrings; j++)
			{
				BIntrinsicForOneEstimation(j, i) = (K(j, i) / 4.0) * pow(dCore(j), 2.0);
				BPluckForOneEstimation(j, i) = ((K(j, i) * 3.0) / 8.0) * pow(deltaHalf(j, i), 2.0);
				f0ForOneEstimation(j, i) = sqrt(T0(j) / mu(j)) / mL0(j, i) / 2.0;
			}
		}

		f0[nombre_estimation] = f0ForOneEstimation;
		BMixed[nombre_estimation] = BIntrinsicForOneEstimation + BPluckForOneEstimation;

	}

	//-----END-LOOP-----//

	MatrixXd BSimulations[NOMBRE_ESTIMATION];
	MatrixXd w0[NOMBRE_ESTIMATION];

	for (int i = 0; i < NOMBRE_ESTIMATION; i++)
	{
		BSimulations[i] = BMixed[i];
		w0[i] = f0[i];
	}

	//------CLASSIFIER-----//

	Model classifier;

	int kk = 0;

	MatrixXd featureMatrixForOneEstimation(2, 78);
	for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; ++nombre_estimation)
	{
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
			meanFeatureMatrix(0, kk) = 0.0;
			meanFeatureMatrix(1, kk) = 0.0;
			for (int nombre_estimation = 0; nombre_estimation < NOMBRE_ESTIMATION; ++nombre_estimation)
			{
				meanFeatureMatrix(0, kk) = meanFeatureMatrix(0, kk)+(featureMatrix[nombre_estimation](0, kk));
				meanFeatureMatrix(1, kk) = meanFeatureMatrix(1, kk)+(featureMatrix[nombre_estimation](1, kk));
			}
			meanFeatureMatrix(0, kk) /= NOMBRE_ESTIMATION;
			meanFeatureMatrix(1, kk) /= NOMBRE_ESTIMATION;
			kk = kk + 1;
		}
	}
	classifier.Mu = meanFeatureMatrix;

	kk = 0;
	VectorXd W(78);
	MatrixXd covMatrixTemp(NOMBRE_ESTIMATION, 2);
	for (int i = 0; i <= numFrets; i++)
	{
		for (int j = numStrings - 1; j >= 0; j--)
		{
			
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
	myfile.close(); // for random
	return classifier;

}
