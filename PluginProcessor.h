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
#include <thread>        
#include <mutex>         
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
#define NBFRET 12
#define NBSTRING 6
#define NOMBRE_ESTIMATION 500

struct Model;

struct result
{
	std::vector<bool> isPlayed;
	std::vector<int> Fret;
	std::vector<int> String;
	std::vector<float> Onset;
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

	void pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size, float onsetValue) noexcept;
	void drawNextFrameOfSpectrum();
	bool getNextFFTBlockReady();
	void setNextFFTBlockReady(bool value);
	void drawFrame(Graphics& g);

	//==============================================================================//
	//==================================Onset=======================================//
	//==============================================================================//

	void setThresholdValue(int value);
	void ProcessWithSpectrogramOnsetDetector(const float* data, int bufSize, int sampleRate);

	//==============================================================================//
	//================================Prediction====================================//
	//==============================================================================//
	
	float            PitchEstimator         (const std::vector<float>& segment);
	std::vector<int> BetaEstimator           (const std::vector<float>& segment, float observedPitch, std::vector<float> pitchReference, float& f0Features, float& betaFeatures);
	int              FretStringPrediction    (const std::vector<float>& segment);
	float           PluckPositionPrediction (const std::vector<std::complex<float>>& hilberOutput, float f0Features, float betaFeatures, int fretPlayed);

	//==============================================================================//
	//============================Processing and test===============================//
	//==============================================================================//

	float getAfficheValue();

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

	std::vector<float> storageSparse;
	std::vector<float> storageActual;
	std::vector<float> storagePast;
	std::vector<float> newSegment;
	std::vector<float> onsetScope;//just for working

	bool detectionDone;
	bool onsetdetected;
	 
	std::vector<bool> timeOfOnset;

	//-------second-step-------//

	//--------//pitch//--------//

	int f0LimitsInf; 
	int f0LimitsSup;           // boundaries for f0 search grid in Hz.
	int nbOfHarmonicsInit;     // number of harmonics for initial harmonic estimate(B = 0).

	//-------third--step-------//

	//---//pitch-candidate//---//




	//-------final-step-------//

	//--------//result//--------//


	//-------final-step-------//

	//--------//result//--------//

	int stringPlayed;
	int fretPlayed;

	//-------final-step-------//

	//--------//plot//--------//
	bool scopeDataIsPlayed[SCOPESIZE];//The scopeData float array of size 512 will contain the points to display on the screen.
	int scopeDataFret     [SCOPESIZE];
	int scopeDataString   [SCOPESIZE];
	float scopeDataOnset  [SCOPESIZE];//just for working
	result DataScope;
	result DataScopeFifo;
	bool nextFFTBlockReady;



	float afficheValue;
	Model classifier;
	Eigen::MatrixXd featureMatrix[NOMBRE_ESTIMATION];

	Eigen::MatrixXd w0; //= squeeze(est_f0);
	Eigen::MatrixXd B;  //= squeeze(BCoeff);

	//sweep//
	int trueString;
	int trueFret;

	//input buffer//
	Eigen::VectorXd X;

	//classifier parameters//
	Eigen::VectorXd mu;
	Eigen::MatrixXd C;
	float P;
	Eigen::VectorXd JMatrix;
	float prediction;
	std::vector<float> pitchReference;

	//Initialize estimator / classifier implementation constants
	float       segmentDuration   = 40e-3;       // segment duration in seconds.
	float       LOpen             = 64.3;        // assumed length of all open strings.
	int highNbOfHarmonics = 25;          // assumed high number of harmonics(M >> 1).Used for inharmonic pitch estimationand for estimation of plucking position on the fretted string.

	float       lengthFFT         = pow(2,19);   // Length of  zero - padded FFT.
	float       betaRes           = 1e-5;        // resolution of search grid for B.



	//result





	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Jacode_iiiAudioProcessor)
};









































//void ProcessWithLowCPUDetector(const float* data, int bufSize, int sampleRate);

	/*void threadPredict(std::vector<float> const& newSegment, int& prediction);
	void Jacode_iiiAudioProcessor::OnsetThread();
	void Jacode_iiiAudioProcessor::PredictionThread();
	void Jacode_iiiAudioProcessor::DisplayThread();
	void Jacode_iiiAudioProcessor::RunThread();*/