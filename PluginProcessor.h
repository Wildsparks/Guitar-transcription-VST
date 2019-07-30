/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream> 
#include <algorithm>    
#include <condition_variable>
#include <numeric>
#include <future>

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\HarmonicSummation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PitchCandidate.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\InharmonicSummation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\AmplitudesEstimation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PluckingPositionEstimatorLSD.h"

#define SCOPESIZE 4096
#define NBFRET 21
#define NBSTRING 6
#define NOMBRE_ESTIMATION 500

struct Model;

struct result
{
	std::vector<bool> isPlayed;
	std::vector<int> Fret;
	std::vector<int> String;
	std::vector<double> Onset;
	std::vector<double> Pluck;
};

//==============================================================================
/**
*/
class Jacode_iiiAudioProcessor  : public AudioProcessor                             //héritage de audioprocessor
{
public:                                                                             //déclaration des méthodes et accesseur
    //==============================================================================
    Jacode_iiiAudioProcessor();                                                     //constructor
    ~Jacode_iiiAudioProcessor();                                                    //décontructor

    //=================================METHODE======================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool   acceptsMidi() const override;
    bool   producesMidi() const override;
    bool   isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int   getNumPrograms() override;
    int   getCurrentProgram() override;
    void  setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void  changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

	//==============================================================================//
	//================================Osciloscope===================================//
	//==============================================================================//

	void pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size, double onsetValue, double pluckValue) noexcept;
	void drawNextFrame();
	bool getNextBlockReady();
	void setNextBlockReady(bool value);
	void drawFrame(Graphics& g);

	//==============================================================================//
	//==================================Onset=======================================//
	//==============================================================================//

	void setThresholdValue(int value);
	void ProcessWithSpectrogramOnsetDetector(const double* data, int bufSize, double sampleRate);
	void ProcessWithLowCPUOnsetDetector(const double* data, int bufSize, double sampleRate);
	//==============================================================================//
	//================================Prediction====================================//
	//==============================================================================//
	
	double           PitchEstimator          (const std::vector<double>& segment);
	std::vector<int> BetaEstimator           (const std::vector<double>& segment, double observedPitch, std::vector<double> pitchReference, double& f0Features, double& betaFeatures);
	void             FretStringPrediction(const std::vector<double>& segment, int& Noteprediction, double& pluckValue);
	double           PluckPositionPrediction (const std::vector<std::complex<double>>& hilberOutput, double f0Features, double betaFeatures, int fretPlayed);

	//==============================================================================//
	//============================ == = editor  = == ===============================//
	//==============================================================================//

	double getAfficheValue();
	void setnfftValue(double value);
	void setbetaresValue(double value);
	void setPluckOnOff(bool value);
	void setResPluck(double value);
	void setDistanceMaxPluck(double value);
	void setNBHarmonics(double value);
	void setNBHarmonicsF0(double value);
	void reloadModel(double val1, double val2, double val3, double val4, double val5,
		double val6, double val7, double val8, double val9, double val10,
		double val11, double val12, double val13, double val14, double val15,
		double val16, double val17, double val18, double val19, double val20,
		double val21, double val22);

	//==============================================================================//
	//================================End method====================================//
	//==============================================================================//


private:                                                                            //déclaration des Attributs


	//==============================================================================//
    //================================ATTRIBUT======================================//
	//==============================================================================//

	//-------first--step-------//

	//--------//onset//--------//

	int counterOnset;
	int thresholdValue;

	std::vector<double> storageSparse;
	std::vector<double> storageActual;
	std::vector<double> storagePast;
	std::vector<double> newSegment;
	std::vector<double> onsetScope;//just for working

	bool detectionDone;
	bool onsetdetected;
	 
	std::vector<bool> timeOfOnset;

	//-------second-step-------//

	//------//prepocess//------//

	double       lengthFFT;   // Length of  zero - padded FFT.
	double       segmentDuration;       // segment duration in seconds.

	//-------third--step-------//

	//--------//pitch//--------//

	int f0LimitsInf;
	int f0LimitsSup;           // boundaries for f0 search grid in Hz.
	int nbOfHarmonicsInit;     // number of harmonics for initial harmonic estimate(B = 0).

	//--------fourth-step-------//

	//------//Candidate//-------//



	//--------fifth-step--------//

	//------//inharmonic//------//

	double betaRes;	            // resolution of search grid for B.
	int    highNbOfHarmonics;

	//--------sixth-step--------//

	//------//Prediction//------//

	std::default_random_engine generator;
	std::normal_distribution<double> distribution;
	Model classifier;
	Eigen::VectorXd X;
	Eigen::MatrixXd featureMatrix[NOMBRE_ESTIMATION];
	Eigen::MatrixXd w0; 
	Eigen::MatrixXd B;  
	Eigen::VectorXd mu;
	Eigen::MatrixXd C;
	double P;
	Eigen::VectorXd JMatrix;
	double prediction;
	std::vector<double> pitchReference;
	int trueString;
	int trueFret;

	//-------seventh-step-------//

	//----//pluck position//----//

	double LOpen;               // assumed length of all open strings.
	double distanceMaxPluck;
	double resPluck;

	//-------final-step-------//

	//--------//result//--------//

	int stringPlayed;
	int fretPlayed;
	bool isPlayed;
	double pluckValue;

	//-------final-step-------//

	//--------//plot//--------//

	bool   scopeDataIsPlayed [SCOPESIZE];//The scopeData float array of size 512 will contain the points to display on the screen.
	int    scopeDataFret     [SCOPESIZE];
	int    scopeDataString   [SCOPESIZE];
	double scopeDataOnset    [SCOPESIZE];
	double scopeDataPluck    [SCOPESIZE];
	result DataScopeFifo;
	bool   nextBlockReady;

	//-----miscellaneous-----//
	bool PluckOnOff;
	double afficheValue;
	//--------//---//--------//

	

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Jacode_iiiAudioProcessor)
};
