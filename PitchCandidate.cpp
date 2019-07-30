#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PitchCandidate.h"
#include <limits>
#include <iostream>
#include <D:/ecler/Documents/Cours/Ingenieur_4A/Stage/Jacode_III/Source/package eigen/Eigen/Dense>

#define NUMFRET 22
#define NUMSTRING 6

void PitchCandidate(std::vector<int>& finalCandidate, double observedPitch, std::vector<double> pitchReference, int nbString)
{
	double minVal  (std::numeric_limits<double>::max());
	std::vector<double>distanceInCents(NUMFRET * NUMSTRING);

	for (int ss = 0; ss < NUMSTRING; ++ss)
	{
		for (int ff = 0; ff < NUMFRET; ++ff)
		{
			distanceInCents[(__int64)ss + NUMSTRING * (__int64)ff] = std::abs(1200.0 * log2(pitchReference[(__int64)ss + NUMSTRING * (__int64)ff] / observedPitch));
		
			// trach the minimum
			if (distanceInCents[(__int64)ss + NUMSTRING * (__int64)ff] < minVal)
			{
				minVal = distanceInCents[(__int64)ss + NUMSTRING * (__int64)ff];
			}
		}
	}
	finalCandidate.clear();
	int point(-2);
	for (int ss = 0; ss < NUMSTRING; ++ss)
	{
		for (int ff = 0; ff < NUMFRET; ++ff)
		{
			if (distanceInCents[(__int64)ss + NUMSTRING * (__int64)ff] < minVal + 50.0)
			{
				/*
				std::cout << "fret : " << ff << " string : " << ss << std::endl;
				if (abs(6.0-(ss + NUMSTRING * ff - point)) < 2.0)
				{
					if (distanceInCents[ss + NUMSTRING * ff] < distanceInCents[point])
					{
						std::cout << "new : " << distanceInCents[ss + NUMSTRING * ff] << " old : " << distanceInCents[point] << std::endl;
						std::cout << "these one win : fret : " << ff << " string : " << ss << std::endl;
						finalCandidate.pop_back();
						finalCandidate.push_back(ss + NUMSTRING * ff);
					}
				}
				else
				{*/
					finalCandidate.push_back((__int64)ss + NUMSTRING * (__int64)ff);
				//}
				//point = ss + NUMSTRING * ff;
			}
		}
	}

	//sort(finalCandidate.begin(), finalCandidate.end());
}