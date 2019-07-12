#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PitchCandidate.h"
#include <limits>
#include <iostream>
#include <D:/ecler/Documents/Cours/Ingenieur_4A/Stage/Jacode_III/Source/package eigen/Eigen/Dense>

void PitchCandidate(std::vector<int>& finalCandidate, double observedPitch, std::vector<double> pitchReference, int nbString)
{
	double minVal  (std::numeric_limits<double>::max());
	int    minIndex(0);

	for (int i = 0; i < pitchReference.size(); ++i)
	{
		if (minVal > abs(log(pitchReference[i]) - log(observedPitch)))
		{
			minVal = abs(log(pitchReference[i]) - log(observedPitch));
			minIndex = i;
		}
	}

	int fretFind  (minIndex / nbString);
	int stringFind(minIndex % nbString);

	Eigen::MatrixXi Candidate(3, 2);
	Candidate(0, 0) = stringFind;
	Candidate(0, 1) = fretFind;
	Candidate(1, 0) = 98;//impossible value
	Candidate(1, 1) = 98;
	Candidate(2, 0) = 98;
	Candidate(2, 1) = 98;

	switch (stringFind)
	{
	case 0:
		if (fretFind < 10 && fretFind > 4)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
		}

		if (fretFind > 9)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind + 2;
			Candidate(2, 1) = fretFind - 10;
		}
		break;
	case 1:
		if (fretFind < 5)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 5;
		}

		if (fretFind > 4 && fretFind < 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind - 1;
			Candidate(2, 1) = fretFind + 5;
		}

		if (fretFind > 7 && fretFind < 10)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
		}

		if (fretFind > 9)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind + 2;
			Candidate(2, 1) = fretFind - 10;
		}

		break;

	case 2:

		if (fretFind < 5)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 5;
		}
		if (fretFind < 3)
		{
			Candidate(2, 0) = stringFind - 2;
			Candidate(2, 1) = fretFind + 10;
		}

		if (fretFind > 4 && fretFind < 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind - 1;
			Candidate(2, 1) = fretFind + 5;
		}

		if (fretFind > 7 && fretFind < 10)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
		}

		if (fretFind > 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind + 2;
			Candidate(2, 1) = fretFind - 9;
		}

		break;

	case 3:

		if (fretFind < 4)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 5;
		}

		if (fretFind < 3)
		{
			Candidate(2, 0) = stringFind - 2;
			Candidate(2, 1) = fretFind + 10;
		}

		if (fretFind > 3 && fretFind < 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 4;
			Candidate(2, 0) = stringFind - 1;
			Candidate(2, 1) = fretFind + 5;
		}

		if (fretFind > 7 && fretFind < 10)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 4;
		}

		if (fretFind > 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 4;
			Candidate(2, 0) = stringFind + 2;
			Candidate(2, 1) = fretFind - 9;
		}

		break;

	case 4:

		if (fretFind < 5)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 4;
		}


		if (fretFind < 4)
		{
			Candidate(2, 0) = stringFind - 2;
			Candidate(2, 1) = fretFind + 9;
		}

		if (fretFind > 4 && fretFind < 9)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
			Candidate(2, 0) = stringFind - 1;
			Candidate(2, 1) = fretFind + 4;
		}

		if (fretFind > 8)
		{
			Candidate(1, 0) = stringFind + 1;
			Candidate(1, 1) = fretFind - 5;
		}

		break;

	case 5:

		if (fretFind < 4)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 5;
			Candidate(2, 0) = stringFind - 2;
			Candidate(2, 1) = fretFind + 9;
		}

		if (fretFind > 3 && fretFind < 8)
		{
			Candidate(1, 0) = stringFind - 1;
			Candidate(1, 1) = fretFind + 5;
		}

		break;

	default:

		break;
	}

	finalCandidate.resize(3);
	if (Candidate(2, 0) == 98)
	{
		finalCandidate.pop_back();
	}
	if (Candidate(1, 0) == 98)
	{
		finalCandidate.pop_back();
	}

	for (int j = 0; j < finalCandidate.size(); ++j)
	{
		finalCandidate[j] = Candidate(j, 0)+Candidate(j, 1)*6;
	}
	//sort(finalCandidate.begin(), finalCandidate.end());
}